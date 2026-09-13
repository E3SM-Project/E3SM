module RtmHorizRemap

!-----------------------------------------------------------------------
! !MODULE: RtmHorizRemap
!
! !DESCRIPTION:
! MOSART wrapper around the shared horizontal remapping infrastructure
! (share/util/shr_horiz_remap_col_mod.F90), used to write a history tape
! directly on a regular lat-lon target grid instead of the native runoff grid.
!
! This is the river-routing half of the FME (Full Model Emulation) online
! output path; EAM, ELM and the MPAS components have equivalent wrappers.
! Writing the remapped tape online removes the offline regrid step from the
! ROF emulation training pipeline.
!
! A tape is remapped by setting rtmhist_horiz_remap_file(t) in user_nl_mosart
! to an ESMF/SCRIP-format mapping file whose source grid is the full
! (rtmlon x rtmlat) runoff grid -- e.g. map_r05_to_gaussian_180by360_shifted.nc.
!
! Note on conservation: the masked remap produces an area-weighted *mean* of
! the source cells covering each target cell. For a flux density (mm/s, the
! units MOSART history fields carry) that is the right operation. Channel
! volumes and discharges (m3/s) are extensive quantities whose sum is not
! preserved by an average, so a remapped MOSART tape is for emulator forcing
! and diagnostics, not for closing a water budget on the target grid.
!
! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_horiz_remap_col_mod, only : shr_horiz_remap_col_t, SHR_FILL_VALUE
  use shr_sys_mod            , only : shr_sys_abort
  use RtmVar                 , only : iulog, rtmlon, rtmlat
  use RtmSpmd                , only : masterproc, mpicom_rof, iam, npes
  use RunoffMod              , only : rtmCTL
  use RtmIO

  implicit none
  save
  private

! !PUBLIC MEMBER FUNCTIONS:
  public :: RtmHorizRemapInit        ! Read the map file and build the comm pattern
  public :: RtmHorizRemapActive      ! Is remapping active for this tape?
  public :: RtmHorizRemapDims        ! Target grid nlon / nlat
  public :: RtmHorizRemapCoords      ! Target grid lon / lat coordinate arrays
  public :: RtmHorizRemapWriteField  ! Remap one field and write it with PIO
  public :: RtmHorizRemapWriteConst  ! Remap and write a time-constant field

! !PUBLIC DATA:
  public :: SHR_FILL_VALUE

! !PRIVATE DATA:
  integer, parameter :: max_remap_tapes = 3   ! must be >= RtmHistFile's max_tapes
  type(shr_horiz_remap_col_t) :: remap(max_remap_tapes)
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

  subroutine RtmHorizRemapInit(t, mapfile)

    ! !DESCRIPTION:
    ! Initialize remapping for history tape t from the given mapping file.

    ! !ARGUMENTS:
    implicit none
    integer         , intent(in) :: t        ! history tape index
    character(len=*), intent(in) :: mapfile  ! mapping file path

    ! !LOCAL VARIABLES:
    integer :: nr, ier
    integer, allocatable :: my_gindex(:)
    character(len=512) :: remap_errmsg
    character(len=*), parameter :: subname = 'RtmHorizRemapInit'
    !-----------------------------------------------------

    if (t < 1 .or. t > max_remap_tapes) then
       call shr_sys_abort(subname//' ERROR: history tape index out of range')
    end if

    if (masterproc) then
       write(iulog,*) trim(subname),': tape ',t,' will be remapped using ',trim(mapfile)
    end if

    ! rtmCTL%gindex(nr) is the index of local cell nr in the global
    ! (rtmlon x rtmlat) grid, which is the numbering the mapping file's source
    ! side uses.
    allocate(my_gindex(rtmCTL%endr - rtmCTL%begr + 1))
    do nr = rtmCTL%begr, rtmCTL%endr
       my_gindex(nr - rtmCTL%begr + 1) = rtmCTL%gindex(nr)
    end do

    call remap(t)%init(trim(mapfile), my_gindex, rtmlon*rtmlat, &
         mpicom_rof, iam, npes, pio_subsystem, ier, remap_errmsg)
    deallocate(my_gindex)

    if (ier /= 0) then
       if (masterproc) write(iulog,*) trim(subname),' ERROR: ',trim(remap_errmsg)
       call shr_sys_abort(subname//' ERROR: '//trim(remap_errmsg))
    end if

    if (masterproc) then
       write(iulog,*) trim(subname),': tape ',t,' target grid is ', &
            remap(t)%nlon(),' x ',remap(t)%nlat(),' (lon x lat)'
    end if

  end subroutine RtmHorizRemapInit

!-----------------------------------------------------------------------

  logical function RtmHorizRemapActive(t)

    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t
    !-----------------------------------------------------

    if (t < 1 .or. t > max_remap_tapes) then
       RtmHorizRemapActive = .false.
    else
       RtmHorizRemapActive = remap(t)%is_active()
    end if

  end function RtmHorizRemapActive

!-----------------------------------------------------------------------

  subroutine RtmHorizRemapDims(t, nlon, nlat)

    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: t
    integer, intent(out) :: nlon, nlat
    !-----------------------------------------------------

    nlon = remap(t)%nlon()
    nlat = remap(t)%nlat()

  end subroutine RtmHorizRemapDims

!-----------------------------------------------------------------------

  subroutine RtmHorizRemapCoords(t, lon, lat)

    ! !DESCRIPTION:
    ! Return the target grid's 1-d coordinate arrays. Allocated here.

    ! !ARGUMENTS:
    implicit none
    integer,               intent(in)  :: t
    real(r8), allocatable, intent(out) :: lon(:)
    real(r8), allocatable, intent(out) :: lat(:)
    !-----------------------------------------------------

    allocate(lon(remap(t)%nlon()))
    allocate(lat(remap(t)%nlat()))
    lon(:) = remap(t)%shared%lon(:)
    lat(:) = remap(t)%shared%lat(:)

  end subroutine RtmHorizRemapCoords

!-----------------------------------------------------------------------

  subroutine RtmHorizRemapWriteField(t, ncid, varname, fld_local, nt, data_type)

    ! !DESCRIPTION:
    ! Remap one history field onto the target grid and write it as time
    ! record nt.

    ! !USES:
    use RtmVar, only : spval
    use pio   , only : pio_setframe, PIO_OFFSET_KIND

    ! !ARGUMENTS:
    implicit none
    integer          , intent(in)    :: t             ! history tape index
    type(file_desc_t), intent(inout) :: ncid          ! open history file
    character(len=*) , intent(in)    :: varname       ! variable name
    real(r8)         , intent(in)    :: fld_local(:)  ! (local runoff cells)
    integer          , intent(in)    :: nt            ! time record index
    integer          , intent(in)    :: data_type     ! PIO type of the variable

    ! !LOCAL VARIABLES:
    real(r8), allocatable :: fld_in(:,:), fld_out(:,:)
    integer               :: varid
    type(var_desc_t)      :: vardesc
    !-----------------------------------------------------

    allocate(fld_in(size(fld_local), 1))
    fld_in(:,1) = fld_local(:)

    call remap(t)%remap_field(fld_in, 1, fld_out, fillval=spval)
    deallocate(fld_in)

    call ncd_inqvid(ncid, varname, varid, vardesc)
    call pio_setframe(ncid, vardesc, int(nt, kind=PIO_OFFSET_KIND))
    call remap(t)%write_field(ncid, vardesc, fld_out, 1, data_type)

    deallocate(fld_out)

  end subroutine RtmHorizRemapWriteField

!-----------------------------------------------------------------------

  subroutine RtmHorizRemapWriteConst(t, ncid, varname, fld_local, data_type, &
       missing_as_zero)

    ! !DESCRIPTION:
    ! Remap and write a time-constant field (no time dimension, so no frame
    ! is set). Pass missing_as_zero for a coverage-like quantity such as the
    ! runoff mask, where cells outside the source decomposition must count as
    ! genuine zeros rather than be excluded from the average.

    ! !ARGUMENTS:
    implicit none
    integer          , intent(in)    :: t
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: varname
    real(r8)         , intent(in)    :: fld_local(:)
    integer          , intent(in)    :: data_type
    logical, optional, intent(in)    :: missing_as_zero

    ! !LOCAL VARIABLES:
    real(r8), allocatable :: fld_in(:,:), fld_out(:,:)
    integer               :: varid
    logical               :: lzero
    type(var_desc_t)      :: vardesc
    !-----------------------------------------------------

    lzero = .false.
    if (present(missing_as_zero)) lzero = missing_as_zero

    allocate(fld_in(size(fld_local), 1))
    fld_in(:,1) = fld_local(:)

    call remap(t)%remap_field(fld_in, 1, fld_out, missing_as_zero=lzero)
    deallocate(fld_in)

    call ncd_inqvid(ncid, varname, varid, vardesc)
    call remap(t)%write_field(ncid, vardesc, fld_out, 1, data_type)

    deallocate(fld_out)

  end subroutine RtmHorizRemapWriteConst

end module RtmHorizRemap
