module restCompactMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Support for compact restart files (restart_file_type = 'compact'), in which
  ! column- and pft-level data are written only for active columns and pfts.
  !
  ! Compact files have the dimensions namec_compact and namep_compact next to the
  ! full column and pft dimensions. Column (pft) variables are defined on the
  ! compact dimension, except the subgrid weights and metadata written by
  ! subgridRestMod, which stay on the full dimensions. The points of the compact
  ! dimension are the active columns (pfts) in the same order as on the full
  ! dimension, so the layout does not depend on the domain decomposition.
  !
  ! This module builds the map from local columns (pfts) to the compact
  ! dimensions that ncdio_pio uses for variables on those dimensions, and
  ! handles the file metadata that describes the compact layout.
  !
  ! !USES:
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use abortutils         , only : endrun
  use spmdMod            , only : masterproc, mpicom, MPI_INTEGER, MPI_SUM
  use elm_varctl         , only : iulog
  use elm_varcon         , only : grlnd, nameg, namec, namep, namec_compact, namep_compact
  use decompMod          , only : bounds_type, BOUNDS_LEVEL_PROC, gsMap_lnd_gdc2glo
  use decompMod          , only : col_compact_gindex, pft_compact_gindex
  use decompMod          , only : numc_compact, nump_compact, compact_map_epoch
  use mct_mod            , only : mct_gsMap_gsize
  use spmdGathScatMod    , only : gather_data_to_master, scatter_data_from_master
  use ncdio_pio          , only : file_desc_t, ncd_int, ncd_global, ncd_defdim, ncd_inqdid, ncd_inqdlen
  use ncdio_pio          , only : ncd_io, ncd_putatt, ncd_set_compact
  use GetGlobalValuesMod , only : GetGlobalIndexArray
  use ColumnType         , only : col_pp
  use VegetationType     , only : veg_pp
  use restUtilMod
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: restCompact_build_map  ! build the compact map from the current active flags
  public :: restCompact_dimset     ! define the compact dimensions on a file being created
  public :: restCompact_ids        ! define/write/check the ids of the points on the compact dimensions
  public :: restCompact_read_map   ! build the compact map from the active flags on a compact file
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: build_level_map
  private :: set_map
  private :: ordinal_in_gridcell

  character(len=*), parameter, public :: restart_file_type_attname = 'restart_file_type'
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine restCompact_build_map(bounds)
    !
    ! !DESCRIPTION:
    ! Build the compact map from the current column and pft active flags. Called
    ! before writing a compact restart file.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer, allocatable :: col_mask(:)
    integer, allocatable :: pft_mask(:)
    !-----------------------------------------------------------------------

    allocate(col_mask(bounds%begc:bounds%endc), pft_mask(bounds%begp:bounds%endp))
    col_mask(:) = 0
    pft_mask(:) = 0
    where (col_pp%active(bounds%begc:bounds%endc)) col_mask = 1
    where (veg_pp%active(bounds%begp:bounds%endp)) pft_mask = 1

    call set_map(bounds, col_mask, pft_mask)

    deallocate(col_mask, pft_mask)

  end subroutine restCompact_build_map

  !-----------------------------------------------------------------------
  subroutine restCompact_read_map(bounds, ncid, compact)
    !
    ! !DESCRIPTION:
    ! Determine whether a restart file is compact. If it is, build the compact
    ! map from the column and pft active flags on the file (which are on the
    ! full dimensions), and check it against the ids on the file.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds   ! proc-level bounds
    type(file_desc_t), intent(inout) :: ncid     ! netcdf id
    logical          , intent(out)   :: compact  ! true => file is compact
    !
    ! !LOCAL VARIABLES:
    integer          :: dimid
    integer          :: dimlen
    logical          :: pft_dim_exists
    logical          :: readvar
    integer, pointer :: col_mask(:)
    integer, pointer :: pft_mask(:)

    character(len=*), parameter :: subname = 'restCompact_read_map'
    !-----------------------------------------------------------------------

    call ncd_inqdid(ncid, namec_compact, dimid, dimexist=compact)
    call ncd_inqdid(ncid, namep_compact, dimid, dimexist=pft_dim_exists)
    if (compact .neqv. pft_dim_exists) then
       call endrun(msg=subname//' ERROR: restart file has only one of the dimensions '// &
            trim(namec_compact)//' and '//trim(namep_compact)//errMsg(__FILE__, __LINE__))
    end if
    if (.not. compact) return

    allocate(col_mask(bounds%begc:bounds%endc), pft_mask(bounds%begp:bounds%endp))

    call ncd_io(ncid=ncid, varname='cols1d_active', flag='read', data=col_mask, &
         dim1name=namec, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=subname//' ERROR: cols1d_active not found on compact restart file'// &
            errMsg(__FILE__, __LINE__))
    end if
    call ncd_io(ncid=ncid, varname='pfts1d_active', flag='read', data=pft_mask, &
         dim1name=namep, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=subname//' ERROR: pfts1d_active not found on compact restart file'// &
            errMsg(__FILE__, __LINE__))
    end if

    call set_map(bounds, col_mask, pft_mask)

    deallocate(col_mask, pft_mask)

    call ncd_inqdid(ncid, namec_compact, dimid)
    call ncd_inqdlen(ncid, dimid, dimlen)
    if (dimlen /= numc_compact) then
       write(iulog,*) subname//' ERROR: ',trim(namec_compact),' = ',dimlen, &
            ' but cols1d_active has ',numc_compact,' active columns'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    call ncd_inqdid(ncid, namep_compact, dimid)
    call ncd_inqdlen(ncid, dimid, dimlen)
    if (dimlen /= nump_compact) then
       write(iulog,*) subname//' ERROR: ',trim(namep_compact),' = ',dimlen, &
            ' but pfts1d_active has ',nump_compact,' active pfts'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    call restCompact_ids(bounds, ncid, flag='read')

    if (masterproc) then
       write(iulog,*) 'Reading compact restart file: ', numc_compact, ' columns and ', &
            nump_compact, ' pfts'
    end if

  end subroutine restCompact_read_map

  !-----------------------------------------------------------------------
  subroutine restCompact_dimset(ncid)
    !
    ! !DESCRIPTION:
    ! Define the compact dimensions on a file being created, and register the
    ! file so that its column and pft variables are defined on them. The compact
    ! map must have been built.
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid     ! netcdf id
    !
    ! !LOCAL VARIABLES:
    integer :: dimid

    character(len=*), parameter :: subname = 'restCompact_dimset'
    !-----------------------------------------------------------------------

    if (compact_map_epoch == 0) then
       call endrun(msg=subname//' ERROR: compact map has not been built'//errMsg(__FILE__, __LINE__))
    end if

    call ncd_defdim(ncid, namec_compact, numc_compact, dimid)
    call ncd_defdim(ncid, namep_compact, nump_compact, dimid)
    call ncd_putatt(ncid, ncd_global, restart_file_type_attname, 'compact')
    call ncd_set_compact(ncid, .true.)

  end subroutine restCompact_dimset

  !-----------------------------------------------------------------------
  subroutine restCompact_ids(bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Ids of the points on the compact dimensions: the global gridcell index and
    ! the 1-based position within that gridcell (in the full column/pft order).
    ! They describe the compact layout for offline tools. On read they are
    ! checked against the map built from the active flags; this also checks that
    ! reading leaves the columns (pfts) not on the file unchanged.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds  ! proc-level bounds
    type(file_desc_t), intent(inout) :: ncid    ! netcdf id
    character(len=*) , intent(in)    :: flag    ! 'define', 'write' or 'read'
    !
    ! !LOCAL VARIABLES:
    integer          :: c, p
    logical          :: readvar
    integer, pointer :: icarr(:)       ! temporary
    integer, pointer :: iparr(:)       ! temporary
    integer, pointer :: col_grc(:)     ! expected gridcell index of each column
    integer, pointer :: col_ord(:)     ! expected position of each column in its gridcell
    integer, pointer :: pft_grc(:)     ! expected gridcell index of each pft
    integer, pointer :: pft_ord(:)     ! expected position of each pft in its gridcell
    integer, parameter :: unset = -999 ! value of points not on the file after read

    character(len=*), parameter :: subname = 'restCompact_ids'
    !-----------------------------------------------------------------------

    allocate(icarr(bounds%begc:bounds%endc), iparr(bounds%begp:bounds%endp))
    allocate(col_grc(bounds%begc:bounds%endc), col_ord(bounds%begc:bounds%endc))
    allocate(pft_grc(bounds%begp:bounds%endp), pft_ord(bounds%begp:bounds%endp))

    col_grc = GetGlobalIndexArray(col_pp%gridcell(bounds%begc:bounds%endc), bounds%begc, bounds%endc, &
         elmlevel=nameg)
    call ordinal_in_gridcell(bounds, col_pp%gridcell(bounds%begc:bounds%endc), col_ord)
    pft_grc = GetGlobalIndexArray(veg_pp%gridcell(bounds%begp:bounds%endp), bounds%begp, bounds%endp, &
         elmlevel=nameg)
    call ordinal_in_gridcell(bounds, veg_pp%gridcell(bounds%begp:bounds%endp), pft_ord)

    if (flag == 'read') then
       icarr(:) = unset
    else
       icarr(:) = col_grc(:)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_active_gridcell_index', xtype=ncd_int, &
         dim1name=namec,                                                                         &
         long_name='gridcell index of column in compact restart file',                          &
         interpinic_flag='skip', readvar=readvar, data=icarr)
    if (flag == 'read') call check_ids('cols1d_active_gridcell_index', icarr, col_grc, col_compact_gindex)

    if (flag == 'read') then
       icarr(:) = unset
    else
       icarr(:) = col_ord(:)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_active_ordinal', xtype=ncd_int, &
         dim1name=namec,                                                                  &
         long_name='1-based position of column within its gridcell',                     &
         interpinic_flag='skip', readvar=readvar, data=icarr)
    if (flag == 'read') call check_ids('cols1d_active_ordinal', icarr, col_ord, col_compact_gindex)

    if (flag == 'read') then
       iparr(:) = unset
    else
       iparr(:) = pft_grc(:)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_active_gridcell_index', xtype=ncd_int, &
         dim1name=namep,                                                                         &
         long_name='gridcell index of pft in compact restart file',                             &
         interpinic_flag='skip', readvar=readvar, data=iparr)
    if (flag == 'read') call check_ids('pfts1d_active_gridcell_index', iparr, pft_grc, pft_compact_gindex)

    if (flag == 'read') then
       iparr(:) = unset
    else
       iparr(:) = pft_ord(:)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_active_ordinal', xtype=ncd_int, &
         dim1name=namep,                                                                  &
         long_name='1-based position of pft within its gridcell',                        &
         interpinic_flag='skip', readvar=readvar, data=iparr)
    if (flag == 'read') call check_ids('pfts1d_active_ordinal', iparr, pft_ord, pft_compact_gindex)

    deallocate(icarr, iparr, col_grc, col_ord, pft_grc, pft_ord)

  contains

    subroutine check_ids(varname, values, expected, gindex)
      ! Points on the file must match the expected ids; the others must be unchanged
      character(len=*), intent(in) :: varname
      integer         , intent(in) :: values(:)
      integer         , intent(in) :: expected(:)
      integer         , intent(in) :: gindex(:)
      integer :: i

      do i = 1, size(values)
         if (gindex(i) > 0) then
            if (values(i) /= expected(i)) then
               write(iulog,*) subname//' ERROR: ',trim(varname),' at compact index ',gindex(i), &
                    ' is ',values(i),', expected ',expected(i)
               call endrun(msg=errMsg(__FILE__, __LINE__))
            end if
         else if (values(i) /= unset) then
            write(iulog,*) subname//' ERROR: reading ',trim(varname), &
                 ' changed a point that is not on the compact file'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if
      end do

    end subroutine check_ids

  end subroutine restCompact_ids

  !-----------------------------------------------------------------------
  subroutine set_map(bounds, col_mask, pft_mask)
    !
    ! !DESCRIPTION:
    ! Build the compact column and pft maps from the given masks.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds                     ! proc-level bounds
    integer          , intent(in) :: col_mask(bounds%begc:)     ! 1 => column is on compact files
    integer          , intent(in) :: pft_mask(bounds%begp:)     ! 1 => pft is on compact files
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'set_map'
    !-----------------------------------------------------------------------

    if (bounds%level /= BOUNDS_LEVEL_PROC) then
       call endrun(msg=subname//' ERROR: expect proc-level bounds'//errMsg(__FILE__, __LINE__))
    end if

    call build_level_map(bounds, col_pp%gridcell(bounds%begc:bounds%endc), &
         col_mask(bounds%begc:bounds%endc), col_compact_gindex, numc_compact)
    call build_level_map(bounds, veg_pp%gridcell(bounds%begp:bounds%endp), &
         pft_mask(bounds%begp:bounds%endp), pft_compact_gindex, nump_compact)

    ! Invalidates io descriptors built for the previous map
    compact_map_epoch = compact_map_epoch + 1

  end subroutine set_map

  !-----------------------------------------------------------------------
  subroutine build_level_map(bounds, gridcell, mask, gindex, gsize)
    !
    ! !DESCRIPTION:
    ! Compact map for one subgrid level (columns or pfts). The compact index of
    ! a point is its position among the masked points in global order: gridcells
    ! in global order, then the points of each gridcell in local order. This is
    ! the order of the full dimension (see decompInit_gtlcp), restricted to the
    ! masked points.
    !
    ! !ARGUMENTS:
    type(bounds_type)   , intent(in)    :: bounds       ! proc-level bounds
    integer             , intent(in)    :: gridcell(:)  ! local gridcell index of each local point
    integer             , intent(in)    :: mask(:)      ! 1 => point is on compact files
    integer, allocatable, intent(inout) :: gindex(:)    ! compact index of each local point; 0 => not on file
    integer             , intent(out)   :: gsize        ! global number of masked points
    !
    ! !LOCAL VARIABLES:
    integer          :: i, g, n, ng
    integer          :: val1, val2
    integer          :: nlocal
    integer          :: ier
    integer, pointer :: gcount(:)     ! number of masked points in each local gridcell
    integer, pointer :: gstart(:)     ! compact index of the first masked point of each local gridcell
    integer, pointer :: arrayglob(:)  ! gridcell values in global order (master only)
    !-----------------------------------------------------------------------

    allocate(gcount(bounds%begg:bounds%endg), gstart(bounds%begg:bounds%endg))
    gcount(:) = 0
    gstart(:) = 0
    do i = 1, size(mask)
       if (mask(i) == 1) then
          g = gridcell(i)
          gcount(g) = gcount(g) + 1
       end if
    end do

    ! Gather the counts to the master in global gridcell order, turn them into
    ! start indices, and scatter these back

    ng = mct_gsMap_gsize(gsMap_lnd_gdc2glo)
    if (masterproc) then
       allocate(arrayglob(ng))
    else
       allocate(arrayglob(1))
    end if
    arrayglob(:) = 0
    call gather_data_to_master(gcount, arrayglob, grlnd)
    if (masterproc) then
       val1 = arrayglob(1)
       arrayglob(1) = 1
       do n = 2,ng
          val2 = arrayglob(n)
          arrayglob(n) = arrayglob(n-1) + val1
          val1 = val2
       enddo
    endif
    call scatter_data_from_master(gstart, arrayglob, grlnd)
    deallocate(arrayglob)

    if (allocated(gindex)) deallocate(gindex)
    allocate(gindex(size(mask)))
    gcount(:) = 0
    do i = 1, size(mask)
       if (mask(i) == 1) then
          g = gridcell(i)
          gindex(i) = gstart(g) + gcount(g)
          gcount(g) = gcount(g) + 1
       else
          gindex(i) = 0
       end if
    end do

    nlocal = count(mask(:) == 1)
    call mpi_allreduce(nlocal, gsize, 1, MPI_INTEGER, MPI_SUM, mpicom, ier)

    deallocate(gcount, gstart)

  end subroutine build_level_map

  !-----------------------------------------------------------------------
  subroutine ordinal_in_gridcell(bounds, gridcell, ordinal)
    !
    ! !DESCRIPTION:
    ! 1-based position of each local point among all points of its gridcell.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)  :: bounds       ! proc-level bounds
    integer          , intent(in)  :: gridcell(:)  ! local gridcell index of each local point
    integer          , intent(out) :: ordinal(:)   ! position within gridcell
    !
    ! !LOCAL VARIABLES:
    integer              :: i, g
    integer, allocatable :: ioff(:)
    !-----------------------------------------------------------------------

    allocate(ioff(bounds%begg:bounds%endg))
    ioff(:) = 0
    do i = 1, size(gridcell)
       g = gridcell(i)
       ioff(g) = ioff(g) + 1
       ordinal(i) = ioff(g)
    end do
    deallocate(ioff)

  end subroutine ordinal_in_gridcell

end module restCompactMod
