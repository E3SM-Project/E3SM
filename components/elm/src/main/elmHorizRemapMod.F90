module elmHorizRemapMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! ELM wrapper around the shared horizontal remapping infrastructure
  ! (share/util/shr_horiz_remap_col_mod.F90), used to write a history tape
  ! directly on a regular lat-lon target grid instead of the native land grid.
  !
  ! This is the ELM half of the FME (Full Model Emulation) online output path;
  ! EAM has an equivalent wrapper in components/eam/src/control/horiz_remap_mod.F90
  ! and MPAS-O / MPAS-SI have theirs in their analysis members. Writing the
  ! remapped tape online removes the offline regrid step from the ELM/ROF
  ! emulation training pipeline.
  !
  ! A tape is remapped by setting hist_horiz_remap_file(t) in user_nl_elm to an
  ! ESMF/SCRIP-format mapping file whose source grid is the full (ni x nj) land
  ! grid -- e.g. map_r05_to_gaussian_180by360_shifted.nc. Only tapes written on
  ! the gridcell decomposition (hist_dov2xy = .true.) can be remapped.
  !
  ! ELM owns only the land cells of that source grid, so the shared layer's
  ! missing-column handling does the work: ocean cells contribute nothing to
  ! the masked average, so a coastal target cell reports the mean over its land
  ! part and a target cell with no land at all gets the fill value.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_horiz_remap_col_mod, only : shr_horiz_remap_col_t, SHR_FILL_VALUE
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use elm_varctl          , only : iulog
  use spmdMod             , only : masterproc, mpicom, iam, npes
  use ncdio_pio           , only : file_desc_t, var_desc_t
  !
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: elm_horiz_remap_init        ! Read the map file and build the comm pattern
  public :: elm_horiz_remap_active      ! Is remapping active for this tape?
  public :: elm_horiz_remap_dims        ! Target grid nlon / nlat
  public :: elm_horiz_remap_coords      ! Target grid lon / lat coordinate arrays
  public :: elm_horiz_remap_write_field ! Remap one field and write it with PIO
  public :: elm_horiz_remap_write_landfrac ! Remap and write the land fraction
  !
  ! !PUBLIC DATA:
  public :: SHR_FILL_VALUE
  !
  ! !PRIVATE DATA:
  integer, parameter :: max_remap_tapes = 6   ! must be >= histFileMod's max_tapes
  type(shr_horiz_remap_col_t) :: remap(max_remap_tapes)
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine elm_horiz_remap_init(t, mapfile)
    !
    ! !DESCRIPTION:
    ! Initialize remapping for history tape t from the given mapping file.
    !
    ! !USES:
    use domainMod  , only : ldomain
    use decompMod  , only : ldecomp, get_proc_bounds
    use ncdio_pio  , only : pio_subsystem
    !
    ! !ARGUMENTS:
    integer         , intent(in) :: t        ! history tape index
    character(len=*), intent(in) :: mapfile  ! mapping file path
    !
    ! !LOCAL VARIABLES:
    integer :: begg, endg, g, ier
    integer, allocatable :: my_gindex(:)
    character(len=512) :: remap_errmsg
    character(len=*), parameter :: subname = 'elm_horiz_remap_init'
    !-----------------------------------------------------------------------

    if (t < 1 .or. t > max_remap_tapes) then
       call endrun(msg=subname//' ERROR: history tape index out of range'// &
            errMsg(__FILE__, __LINE__))
    end if

    if (masterproc) then
       write(iulog,*) subname,': tape ',t,' will be remapped using ',trim(mapfile)
    end if

    call get_proc_bounds(begg, endg)

    ! ldecomp%gdc2glo(g) is the index of gridcell g in the global (ni x nj)
    ! grid, which is the numbering the mapping file's source side uses.
    allocate(my_gindex(endg-begg+1))
    do g = begg, endg
       my_gindex(g-begg+1) = ldecomp%gdc2glo(g)
    end do

    call remap(t)%init(trim(mapfile), my_gindex, ldomain%ns, &
         mpicom, iam, npes, pio_subsystem, ier, remap_errmsg)
    deallocate(my_gindex)

    if (ier /= 0) then
       if (masterproc) write(iulog,*) subname,' ERROR: ',trim(remap_errmsg)
       call endrun(msg=subname//' ERROR: '//trim(remap_errmsg)//' '// &
            errMsg(__FILE__, __LINE__))
    end if

    if (masterproc) then
       write(iulog,*) subname,': tape ',t,' target grid is ', &
            remap(t)%nlon(),' x ',remap(t)%nlat(),' (lon x lat)'
    end if

  end subroutine elm_horiz_remap_init

  !-----------------------------------------------------------------------
  logical function elm_horiz_remap_active(t)
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t
    !-----------------------------------------------------------------------

    if (t < 1 .or. t > max_remap_tapes) then
       elm_horiz_remap_active = .false.
    else
       elm_horiz_remap_active = remap(t)%is_active()
    end if

  end function elm_horiz_remap_active

  !-----------------------------------------------------------------------
  subroutine elm_horiz_remap_dims(t, nlon, nlat)
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: t
    integer, intent(out) :: nlon, nlat
    !-----------------------------------------------------------------------

    nlon = remap(t)%nlon()
    nlat = remap(t)%nlat()

  end subroutine elm_horiz_remap_dims

  !-----------------------------------------------------------------------
  subroutine elm_horiz_remap_coords(t, lon, lat)
    !
    ! !DESCRIPTION:
    ! Return the target grid's 1-d coordinate arrays. Allocated here.
    !
    ! !ARGUMENTS:
    integer,               intent(in)  :: t
    real(r8), allocatable, intent(out) :: lon(:)
    real(r8), allocatable, intent(out) :: lat(:)
    !-----------------------------------------------------------------------

    allocate(lon(remap(t)%nlon()))
    allocate(lat(remap(t)%nlat()))
    lon(:) = remap(t)%shared%lon(:)
    lat(:) = remap(t)%shared%lat(:)

  end subroutine elm_horiz_remap_coords

  !-----------------------------------------------------------------------
  subroutine elm_horiz_remap_write_field(t, ncid, varname, fld_local, numlev, &
       nt, data_type)
    !
    ! !DESCRIPTION:
    ! Remap fld_local (the gridcell-decomposed history buffer for one field)
    ! onto the target grid and write it as time record nt.
    !
    ! !USES:
    use elm_varcon , only : spval
    use ncdio_pio  , only : ncd_inqvid
    use pio        , only : pio_setframe, PIO_OFFSET_KIND
    !
    ! !ARGUMENTS:
    integer               , intent(in)    :: t          ! history tape index
    type(file_desc_t)     , intent(inout) :: ncid       ! open history file
    character(len=*)      , intent(in)    :: varname    ! variable name
    real(r8)              , intent(in)    :: fld_local(:,:) ! (local gridcells, numlev)
    integer               , intent(in)    :: numlev     ! number of vertical levels
    integer               , intent(in)    :: nt         ! time record index
    integer               , intent(in)    :: data_type  ! PIO type of the variable
    !
    ! !LOCAL VARIABLES:
    real(r8), allocatable :: fld_out(:,:)
    integer               :: varid
    type(var_desc_t)      :: vardesc
    !-----------------------------------------------------------------------

    call remap(t)%remap_field(fld_local, numlev, fld_out, fillval=spval)

    call ncd_inqvid(ncid, varname, varid, vardesc)
    call pio_setframe(ncid, vardesc, int(nt, PIO_OFFSET_KIND))
    call remap(t)%write_field(ncid, vardesc, fld_out, numlev, data_type)

    deallocate(fld_out)

  end subroutine elm_horiz_remap_write_field

  !-----------------------------------------------------------------------
  subroutine elm_horiz_remap_write_landfrac(t, ncid, varname, frac_local, data_type)
    !
    ! !DESCRIPTION:
    ! Remap and write the land fraction. Unlike a prognostic field, the cells
    ! ELM does not own must count as genuine zeros here -- otherwise every
    ! target cell containing any land at all would report a fraction near 1.
    !
    ! !USES:
    use ncdio_pio, only : ncd_inqvid
    !
    ! !ARGUMENTS:
    integer          , intent(in)    :: t
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: varname
    real(r8)         , intent(in)    :: frac_local(:)   ! (local gridcells)
    integer          , intent(in)    :: data_type
    !
    ! !LOCAL VARIABLES:
    real(r8), allocatable :: fld_in(:,:), fld_out(:,:)
    integer               :: varid
    type(var_desc_t)      :: vardesc
    !-----------------------------------------------------------------------

    allocate(fld_in(size(frac_local), 1))
    fld_in(:,1) = frac_local(:)

    call remap(t)%remap_field(fld_in, 1, fld_out, missing_as_zero=.true.)
    deallocate(fld_in)

    call ncd_inqvid(ncid, varname, varid, vardesc)
    call remap(t)%write_field(ncid, vardesc, fld_out, 1, data_type)

    deallocate(fld_out)

  end subroutine elm_horiz_remap_write_landfrac

end module elmHorizRemapMod
