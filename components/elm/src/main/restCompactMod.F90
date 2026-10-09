module restCompactMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Support for compact restart files (restart_file_type = 'compact'), in which
  ! column- and pft-level data are written only for the columns and pfts that
  ! have been active at some time since the base file was read.
  !
  ! The state of a column (pft) is only updated while it is active (except for
  ! the surface albedo variables, see below). A point that has not been active
  ! since the model state was read from a default-layout file (the base file:
  ! finidat, or a default restart file) therefore still has its values from that
  ! file, or its cold start values if there was no such file. Compact files
  ! record the base file, and on read the column/pft variables are first read
  ! from the base file and then from the compact file, so that all points get
  ! the values they would have after reading a default restart file.
  !
  ! Compact files have the dimensions namec_compact and namep_compact next to the
  ! full column and pft dimensions. Column (pft) variables are defined on the
  ! compact dimension, except the subgrid weights and metadata written by
  ! subgridRestMod and the surface albedo variables (computed for inactive points
  ! too), which stay on the full dimensions. The points of the compact
  ! dimension are in the same order as on the full dimension, so the layout
  ! does not depend on the domain decomposition; they are flagged in
  ! cols1d_compact and pfts1d_compact on the full dimensions.
  !
  ! This module tracks which points have been active, builds the map from local
  ! columns (pfts) to the compact dimensions that ncdio_pio uses for variables
  ! on those dimensions, and handles the file metadata that describes the
  ! compact layout and its base file.
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
  use ncdio_pio          , only : ncd_io, ncd_putatt, ncd_set_compact, ncd_compact_suspend
  use ncdio_pio          , only : check_att, ncd_getatt, ncd_pio_openfile, ncd_pio_closefile, ncd_nowrite
  use ncdio_pio          , only : ncd_compact_base_file, ncd_compact_base_set, ncd_compact_base_clear
  use pio                , only : pio_inq_att, PIO_OFFSET_KIND
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
  public :: restCompact_update_ever_active ! add the currently active points to those ever active
  public :: restCompact_build_map  ! build the compact map from the points ever active
  public :: restCompact_dimset     ! define the compact dimensions on a file being created
  public :: restCompact_ids        ! define/write/check the ids of the points on the compact dimensions
  public :: restCompact_read_map   ! determine the layout of a file being read, and set up reading it
  public :: restCompact_read_done  ! close the base file once a compact file has been read
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: build_level_map
  private :: set_map
  private :: ordinal_in_gridcell
  private :: get_global_att

  character(len=*), parameter, public :: restart_file_type_attname = 'restart_file_type'

  ! Global attributes of compact files that identify their base file
  character(len=*), parameter :: base_file_attname    = 'compact_base_file'
  character(len=*), parameter :: base_case_id_attname = 'compact_base_case_id'
  character(len=*), parameter :: base_history_attname = 'compact_base_history'
  character(len=*), parameter :: no_base = 'none'  ! base file attribute for cold start

  ! Points that have been active since the base file was read
  logical, allocatable :: col_ever_active(:)
  logical, allocatable :: pft_ever_active(:)

  ! Base file of compact files written by this run: the last default-layout
  ! restart or initial file read, or no_base for cold start. The case_id and
  ! history (creation time) attributes of the base file are recorded to check
  ! that the same file is read back.
  character(len=512) :: base_file    = no_base
  character(len=256) :: base_case_id = ' '
  character(len=256) :: base_history = ' '
  logical            :: base_open    = .false.  ! true => ncd_compact_base_file is open
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine restCompact_update_ever_active(bounds)
    !
    ! !DESCRIPTION:
    ! Add the currently active columns and pfts to those that have been active
    ! since the base file was read. Must be called whenever the active flags
    ! change: after reading a restart file and after the dynamic subgrid update
    ! of each time step.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !-----------------------------------------------------------------------

    if (.not. allocated(col_ever_active)) then
       allocate(col_ever_active(bounds%begc:bounds%endc), pft_ever_active(bounds%begp:bounds%endp))
       col_ever_active(:) = .false.
       pft_ever_active(:) = .false.
    end if

    col_ever_active(bounds%begc:bounds%endc) = col_ever_active(bounds%begc:bounds%endc) .or. &
         col_pp%active(bounds%begc:bounds%endc)
    pft_ever_active(bounds%begp:bounds%endp) = pft_ever_active(bounds%begp:bounds%endp) .or. &
         veg_pp%active(bounds%begp:bounds%endp)

  end subroutine restCompact_update_ever_active

  !-----------------------------------------------------------------------
  subroutine restCompact_build_map(bounds)
    !
    ! !DESCRIPTION:
    ! Build the compact map from the columns and pfts that have been active since
    ! the base file was read. Called before writing a compact restart file.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer, allocatable :: col_mask(:)
    integer, allocatable :: pft_mask(:)
    !-----------------------------------------------------------------------

    call restCompact_update_ever_active(bounds)

    allocate(col_mask(bounds%begc:bounds%endc), pft_mask(bounds%begp:bounds%endp))
    col_mask(:) = 0
    pft_mask(:) = 0
    where (col_ever_active(bounds%begc:bounds%endc)) col_mask = 1
    where (pft_ever_active(bounds%begp:bounds%endp)) pft_mask = 1

    call set_map(bounds, col_mask, pft_mask)

    deallocate(col_mask, pft_mask)

    if (masterproc) then
       write(iulog,*) 'Compact restart file: ', numc_compact, ' columns and ', &
            nump_compact, ' pfts; base file: ', trim(base_file)
    end if

  end subroutine restCompact_build_map

  !-----------------------------------------------------------------------
  subroutine restCompact_read_map(bounds, ncid, file, compact)
    !
    ! !DESCRIPTION:
    ! Determine whether a restart file is compact, and set up reading it.
    !
    ! A default-layout file becomes the base file of the compact files written
    ! later in the run.
    !
    ! For a compact file, build the compact map from the flags on the file (which
    ! are on the full dimensions), check it against the ids on the file, and
    ! open the base file of the compact file, which is read first for the
    ! column/pft variables on the compact dimensions. Points not on the compact
    ! file are left unchanged if it has no base file (cold start, or a file
    ! written before base files were recorded).
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds   ! proc-level bounds
    type(file_desc_t), intent(inout) :: ncid     ! netcdf id
    character(len=*) , intent(in)    :: file     ! name of the file being read
    logical          , intent(out)   :: compact  ! true => file is compact
    !
    ! !LOCAL VARIABLES:
    integer            :: dimid
    integer            :: dimlen
    logical            :: pft_dim_exists
    logical            :: readvar
    logical            :: exists
    integer, pointer   :: col_mask(:)
    integer, pointer   :: pft_mask(:)
    character(len=512) :: file_base          ! base file recorded on the compact file
    character(len=256) :: file_base_case_id  ! case_id of the base file recorded on the compact file
    character(len=256) :: file_base_history  ! history of the base file recorded on the compact file
    character(len=256) :: case_id            ! case_id attribute of the base file
    character(len=256) :: history            ! history attribute of the base file

    character(len=*), parameter :: subname = 'restCompact_read_map'
    !-----------------------------------------------------------------------

    call ncd_inqdid(ncid, namec_compact, dimid, dimexist=compact)
    call ncd_inqdid(ncid, namep_compact, dimid, dimexist=pft_dim_exists)
    if (compact .neqv. pft_dim_exists) then
       call endrun(msg=subname//' ERROR: restart file has only one of the dimensions '// &
            trim(namec_compact)//' and '//trim(namep_compact)//errMsg(__FILE__, __LINE__))
    end if

    if (allocated(col_ever_active)) deallocate(col_ever_active, pft_ever_active)

    if (.not. compact) then
       ! All points are on this file, so it is the base file of the compact files
       ! written from now on; restFile_read then flags the active points
       base_file = file
       call get_global_att(ncid, 'case_id', base_case_id)
       call get_global_att(ncid, 'history', base_history)
       return
    end if

    allocate(col_mask(bounds%begc:bounds%endc), pft_mask(bounds%begp:bounds%endp))

    ! The points on the file; files written before the base file was recorded
    ! only have the active points
    call ncd_io(ncid=ncid, varname='cols1d_compact', flag='read', data=col_mask, &
         dim1name=namec, readvar=readvar)
    if (.not. readvar) then
       call ncd_io(ncid=ncid, varname='cols1d_active', flag='read', data=col_mask, &
            dim1name=namec, readvar=readvar)
    end if
    if (.not. readvar) then
       call endrun(msg=subname//' ERROR: neither cols1d_compact nor cols1d_active found on compact restart file'// &
            errMsg(__FILE__, __LINE__))
    end if
    call ncd_io(ncid=ncid, varname='pfts1d_compact', flag='read', data=pft_mask, &
         dim1name=namep, readvar=readvar)
    if (.not. readvar) then
       call ncd_io(ncid=ncid, varname='pfts1d_active', flag='read', data=pft_mask, &
            dim1name=namep, readvar=readvar)
    end if
    if (.not. readvar) then
       call endrun(msg=subname//' ERROR: neither pfts1d_compact nor pfts1d_active found on compact restart file'// &
            errMsg(__FILE__, __LINE__))
    end if

    call set_map(bounds, col_mask, pft_mask)

    ! The points on the file are those that had been active since its base file
    ! was read
    allocate(col_ever_active(bounds%begc:bounds%endc), pft_ever_active(bounds%begp:bounds%endp))
    col_ever_active(:) = (col_mask(:) == 1)
    pft_ever_active(:) = (pft_mask(:) == 1)

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

    ! Base file

    call get_global_att(ncid, base_file_attname, file_base)
    if (file_base == ' ') then
       if (masterproc) then
          write(iulog,*) subname//' WARNING: compact restart file has no base file; columns and pfts'
          write(iulog,*) '   not on the file keep their cold start values'
       end if
       file_base = no_base
    end if
    call get_global_att(ncid, base_case_id_attname, file_base_case_id)
    call get_global_att(ncid, base_history_attname, file_base_history)

    if (trim(file_base) /= no_base) then
       inquire(file=trim(file_base), exist=exists)
       if (.not. exists) then
          write(iulog,*) subname//' ERROR: base file of compact restart file not found: ', trim(file_base)
          write(iulog,*) '   The base file holds the state of the columns and pfts not on the'
          write(iulog,*) '   compact file, and must be kept as long as compact files that use it'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       call ncd_pio_openfile(ncd_compact_base_file, trim(file_base), ncd_nowrite)
       base_open = .true.
       call get_global_att(ncd_compact_base_file, 'case_id', case_id)
       call get_global_att(ncd_compact_base_file, 'history', history)
       if (case_id /= file_base_case_id .or. history /= file_base_history) then
          write(iulog,*) subname//' ERROR: base file ', trim(file_base), &
               ' is not the file used for the compact restart file'
          write(iulog,*) '   case_id: expected ', trim(file_base_case_id), ', found ', trim(case_id)
          write(iulog,*) '   history: expected ', trim(file_base_history), ', found ', trim(history)
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       call ncd_inqdid(ncd_compact_base_file, namec_compact, dimid, dimexist=exists)
       if (exists) then
          call endrun(msg=subname//' ERROR: base file '//trim(file_base)// &
               ' is itself a compact restart file'//errMsg(__FILE__, __LINE__))
       end if
       call ncd_compact_base_set(ncid)
       if (masterproc) then
          write(iulog,*) 'Reading columns and pfts not on the compact file from base file ', trim(file_base)
       end if
    end if

    ! Compact files written from now on have the same base file
    base_file    = file_base
    base_case_id = file_base_case_id
    base_history = file_base_history

  end subroutine restCompact_read_map

  !-----------------------------------------------------------------------
  subroutine restCompact_read_done()
    !
    ! !DESCRIPTION:
    ! Close the base file once a compact restart file has been read.
    !-----------------------------------------------------------------------

    if (base_open) then
       call ncd_compact_base_clear()
       call ncd_pio_closefile(ncd_compact_base_file)
       base_open = .false.
    end if

  end subroutine restCompact_read_done

  !-----------------------------------------------------------------------
  subroutine get_global_att(ncid, attname, value)
    !
    ! !DESCRIPTION:
    ! Get a character global attribute; blank if it is not on the file.
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid     ! netcdf id
    character(len=*) , intent(in)    :: attname  ! attribute name
    character(len=*) , intent(out)   :: value    ! attribute value
    !
    ! !LOCAL VARIABLES:
    logical :: found
    integer :: i
    integer :: status
    integer :: att_type
    integer(PIO_OFFSET_KIND) :: att_len

    character(len=*), parameter :: subname = 'get_global_att'
    !-----------------------------------------------------------------------

    value = ' '
    call check_att(ncid, ncd_global, attname, found)
    if (found) then
       ! PIO does not check the length of the value
       status = pio_inq_att(ncid, ncd_global, trim(attname), att_type, att_len)
       if (att_len > len(value)) then
          write(iulog,*) subname//' ERROR: global attribute ',trim(attname),' is longer than ',len(value)
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       call ncd_getatt(ncid, ncd_global, attname, value)
       ! The value is not blank padded beyond the attribute length
       i = index(value, achar(0))
       if (i > 0) value(i:) = ' '
    end if

  end subroutine get_global_att

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
    call ncd_putatt(ncid, ncd_global, base_file_attname, trim(base_file))
    if (trim(base_file) /= no_base) then
       call ncd_putatt(ncid, ncd_global, base_case_id_attname, trim(base_case_id))
       call ncd_putatt(ncid, ncd_global, base_history_attname, trim(base_history))
    end if
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

    ! Flags of the points on the compact dimensions, on the full dimensions. On
    ! read they are read by restCompact_read_map to build the map.

    if (flag /= 'read') then
       call ncd_compact_suspend(.true.)
       icarr(:) = merge(1, 0, col_compact_gindex(:) > 0)
       call restartvar(ncid=ncid, flag=flag, varname='cols1d_compact', xtype=ncd_int, &
            dim1name=namec,                                                          &
            long_name='column is on the compact dimension (1=yes, 0=no)',           &
            interpinic_flag='skip', readvar=readvar, data=icarr)
       iparr(:) = merge(1, 0, pft_compact_gindex(:) > 0)
       call restartvar(ncid=ncid, flag=flag, varname='pfts1d_compact', xtype=ncd_int, &
            dim1name=namep,                                                          &
            long_name='pft is on the compact dimension (1=yes, 0=no)',              &
            interpinic_flag='skip', readvar=readvar, data=iparr)
       call ncd_compact_suspend(.false.)
    end if

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
