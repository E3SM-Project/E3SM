module dshr_dfield_mod

  use ESMF             , only : ESMF_State, ESMF_SUCCESS
  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_sys_mod      , only : shr_sys_abort
  use shr_strdata_mod  , only : shr_strdata_type
  use dshr_methods_mod , only : chkerr, state_getfldptr
  use mct_mod          , only : mct_avect_lsize, mct_avect_indexra

  implicit none
  private

  public :: dshr_dfield_add
  public :: dshr_dfield_copy

  interface dshr_dfield_add
     module procedure dshr_dfield_add_1d
     module procedure dshr_dfield_add_2d
     module procedure dshr_dfield_add_strmfld
  end interface dshr_dfield_add

  ! Linked list node
  type, public :: dfield_type
     ! state data
     real(r8), pointer          :: state_data1d(:) => null()
     real(r8), pointer          :: state_data2d(:,:) => null()
     ! stream data input (always assumed to be 1d for now)
     real(r8), pointer          :: stream_data1d(:) => null()
     ! stream data pointers for 1d state data
     integer                    :: sdat_stream_index = 0
     integer                    :: sdat_avect_index = 0
     ! stream data pointers for 2d state data
     integer, pointer           :: sdat_stream_indices(:) => null()
     integer, pointer           :: sdat_avect_indices(:) => null()
     integer                    :: stream_nflds = 0 ! number of stream field names
     type(dfield_type), pointer :: next => null()
  end type dfield_type

  integer :: iunset = -999
  character(*), parameter :: modName =  "(dshr_dfield_mod)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dshr_dfield_add_strmfld(dfields, sdat, strm_fld, strm_ptr, logunit, masterproc)

    ! input/output variables
    type(dfield_type)      , pointer       :: dfields
    type(shr_strdata_type) , intent(in)    :: sdat
    character(len=*)       , intent(in)    :: strm_fld
    real(r8)               , pointer       :: strm_ptr(:)
    integer, optional      , intent(in)    :: logunit
    logical, optional      , intent(in)    :: masterproc

    ! local variables
    type(dfield_type), pointer :: dfield_new
    integer :: ns, kf
    integer :: lsize, num
    character(len=*), parameter :: subname='(dfield_add_strmfld)'
    ! ----------------------------------------------

    allocate(dfield_new)
    dfield_new%next => dfields
    dfields => dfield_new

    ! determine local size
    lsize = mct_avect_lsize(sdat%avs(1))

    ! initialize dfield_new values
    dfield_new%sdat_stream_index = iunset
    dfield_new%sdat_avect_index = iunset
    dfield_new%stream_data1d => null()

    ! loop over all input streams and ! determine if the strm_fld is in the attribute vector of stream ns
    do ns = 1, sdat%nstreams
       ! determine if the strm_fld is in the attribute vector of stream ns
       kf  = mct_aVect_indexRA(sdat%avs(ns), trim(strm_fld), perrWith='quiet')
       if (kf > 0) then
          dfield_new%sdat_stream_index = ns
          dfield_new%sdat_avect_index = kf
          allocate(dfield_new%stream_data1d(lsize))
          strm_ptr => dfield_new%stream_data1d
          if (present(logunit) .and. present(masterproc)) then
             if (masterproc) then
                write(logunit,*)'(dshr_addfield_add) allocating memory for stream field strm_'//trim(strm_fld)
             end if
          end if
          exit
       end if
    end do

  end subroutine dshr_dfield_add_strmfld

!===============================================================================

  subroutine dshr_dfield_add_1d(dfields, sdat, state_fld, strm_fld, state, &
       state_ptr, strm_ptr, logunit, masterproc, rc)

    ! Set 1d dfield values

    type(dfield_type)      , pointer       :: dfields
    type(shr_strdata_type) , intent(in)    :: sdat
    character(len=*)       , intent(in)    :: state_fld
    character(len=*)       , intent(in)    :: strm_fld
    type(ESMF_State)       , intent(inout) :: state
    real(r8), optional     , pointer       :: state_ptr(:)
    real(r8), optional     , pointer       :: strm_ptr(:)
    integer , optional     , intent(in)    :: logunit
    logical , optional     , intent(in)    :: masterproc
    integer                , intent(out)   :: rc

    ! local variables
    type(dfield_type), pointer :: dfield_new
    integer :: ns, kf, lsize
    character(len=*), parameter :: subname='(dfield_add_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    allocate(dfield_new)
    dfield_new%next => dfields
    dfields => dfield_new

    ! determine local size
    lsize = mct_avect_lsize(sdat%avs(1))

    ! Initialize strm_ptr and state_ptr if it is present
    ! These will be set to valid values if the relevant fields are found
    if (present(strm_ptr)) strm_ptr => null()
    if (present(state_ptr)) state_ptr => null()
    dfield_new%sdat_stream_index = iunset
    dfield_new%sdat_avect_index = iunset

    ! always allocate memory for stream_data and initialize it to 0
    ! note that if the attribute vector is not found in the streams

    ! loop over all input streams
    do ns = 1, sdat%nstreams
       ! determine if the strm_fld is in the attribute vector of stream ns
       kf  = mct_aVect_indexRA(sdat%avs(ns), trim(strm_fld), perrWith='quiet')
       if (kf > 0) then
          dfield_new%sdat_stream_index = ns
          dfield_new%sdat_avect_index = kf
          allocate(dfield_new%stream_data1d(lsize))
          dfield_new%stream_data1d(:) = 0._r8
          if (present(strm_ptr)) then
             strm_ptr => dfield_new%stream_data1d
          end if
          if (present(logunit) .and. present(masterproc)) then
             if (masterproc) then
                write(logunit,*)'(dshr_addfield_add) allocating memory for field '//trim(strm_fld)
             end if
          end if
          exit
       end if
    end do

    ! Set export state array pointer
    call state_getfldptr(State, fldname=trim(state_fld), fldptr1=dfield_new%state_data1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dfield_new%state_data1d = 0.0_r8
    if (present(logunit) .and. present(masterproc)) then
       if (masterproc) then
          write(logunit,*)'(dshr_addfield_add) setting pointer for export state '//trim(state_fld)
       end if
    end if

    ! Return array pointer if argument is present
    if (present(state_ptr)) then
       state_ptr => dfield_new%state_data1d
    end if

  end subroutine dshr_dfield_add_1d

  !===============================================================================

  subroutine dshr_dfield_add_2d(dfields, sdat, state_fld, strm_flds, state, &
       state_ptr, logunit, masterproc, rc)

    ! input/output variables
    type(dfield_type)      , pointer       :: dfields
    type(shr_strdata_type) , intent(in)    :: sdat
    character(len=*)       , intent(in)    :: state_fld
    character(len=*)       , intent(in)    :: strm_flds(:)
    type(ESMF_State)       , intent(inout) :: state
    real(r8), optional     , pointer       :: state_ptr(:,:)
    integer , optional     , intent(in)    :: logunit
    logical , optional     , intent(in)    :: masterproc
    integer                , intent(out)   :: rc

    ! local variables
    type(dfield_type), pointer :: dfield_new
    integer :: n, i, kf, ns, nf
    integer :: nflds, lsize, num
    character(len=*), parameter :: subname='(dfield_add_2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    allocate(dfield_new)
    dfield_new%next => dfields
    dfields => dfield_new

    ! determine stream fldnames array
    nflds = size(strm_flds)
    dfield_new%stream_nflds = nflds
    allocate(dfield_new%sdat_stream_indices(nflds))
    allocate(dfield_new%sdat_avect_indices(nflds))

    ! determine local size
    lsize = mct_avect_lsize(sdat%avs(1))

    ! loop through input array of stream field names
    do nf = 1, nflds
       ! loop through input streams
       do ns = 1, sdat%nstreams
          ! determine if the strm_flds(nf) is in the attribute vector of stream ns
          kf = mct_aVect_indexRA(sdat%avs(ns), trim(strm_flds(nf)), perrWith='quiet')
          if (kf > 0) then
             dfield_new%sdat_stream_indices(nf) = ns
             dfield_new%sdat_avect_indices(nf) = kf
          end if
       end do
    end do

    ! Set export state array pointer
    call state_getfldptr(State, fldname=trim(state_fld), fldptr2=dfield_new%state_data2d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dfield_new%state_data2d(nflds,lsize) = 0._r8

    ! Return array pointer if argument is present
    if (present(state_ptr)) then
       state_ptr => dfield_new%state_data2d
    end if

  end subroutine dshr_dfield_add_2d

  !===============================================================================

  subroutine dshr_dfield_copy(dfields, sdat, rc)

    ! Copy stream data into dfield data type for each element of dfields

    ! input/output variables
    type(dfield_type)      , pointer     :: dfields
    type(shr_strdata_type) , intent(in)  :: sdat
    integer                , intent(out) :: rc

    ! local variables
    type(dfield_type), pointer :: dfield
    integer :: n, i, k
    integer :: avect_index
    integer :: stream_index
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Loop over all dfield entries and fill in stream_data and state_data1d or state_data2d arrays
    dfield => dfields ! note that dfields is the head of the linked list
    do while (associated(dfield))

       ! Map the stream data to the state data
       if (associated(dfield%state_data1d)) then
          stream_index = dfield%sdat_stream_index
          avect_index =  dfield%sdat_avect_index
          if (stream_index /= iunset .and. avect_index /= iunset) then
             dfield%stream_data1d(:) = sdat%avs(stream_index)%rattr(avect_index,:)
             dfield%state_data1d(:) = dfield%stream_data1d(:)
          end if
       else if (associated(dfield%state_data2d)) then
          do i = 1,dfield%stream_nflds
             stream_index = dfield%sdat_stream_indices(i)
             avect_index = dfield%sdat_avect_indices(i)
             dfield%state_data2d(i,:) = sdat%avs(stream_index)%rattr(avect_index,:)
          end do
       else if (associated(dfield%stream_data1d)) then
          stream_index = dfield%sdat_stream_index
          avect_index =  dfield%sdat_avect_index
          if (stream_index /= iunset .and. avect_index /= iunset) then
             dfield%stream_data1d(:) = sdat%avs(stream_index)%rattr(avect_index,:)
          end if
       end if
       dfield => dfield%next

    end do

  end subroutine dshr_dfield_copy

end module dshr_dfield_mod
