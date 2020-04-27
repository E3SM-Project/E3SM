module dshr_dfield_mod

  use ESMF             
  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_sys_mod      , only : shr_sys_abort
  use dshr_strdata_mod , only : shr_strdata_type
  use dshr_methods_mod , only : dshr_state_getfldptr, dshr_fldbun_getfieldn, dshr_field_getfldptr
  use dshr_methods_mod , only : chkerr 

  implicit none
  private

  public :: dshr_dfield_add
  public :: dshr_dfield_copy

  interface dshr_dfield_add
     module procedure dshr_dfield_add_1d
     module procedure dshr_dfield_add_2d
     module procedure dshr_dfield_add_strmfld
  end interface dshr_dfield_add

  ! Note that whereas the data model export state field bundle might have fields
  ! with undistributed dimensions - the stream field bundles only have fields
  ! with no undistributed dimensions

  ! Linked list node
  type, public :: dfield_type
     ! state data
     real(r8), pointer          :: state_data1d(:) => null()
     real(r8), pointer          :: state_data2d(:,:) => null()
     ! stream data input (always assumed to be 1d for now)
     real(r8), pointer          :: stream_data1d(:) => null()
     ! stream data pointers for 1d export state data
     integer                    :: sdat_stream_index = 0
     integer                    :: sdat_fldbun_index = 0
     ! stream data pointers for 2d export state data
     integer, pointer           :: sdat_stream_indices(:) => null()
     integer, pointer           :: sdat_fldbun_indices(:) => null()
     integer                    :: stream_nflds = 0 ! number of stream field names
     character(CS), allocatable :: stream_fldnames(:)
     ! linked list pointer
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
    
    ! Add a dfields element with just stream data info

    ! input/output variables
    type(dfield_type)      , pointer       :: dfields
    type(shr_strdata_type) , intent(in)    :: sdat
    character(len=*)       , intent(in)    :: strm_fld
    real(r8)               , pointer       :: strm_ptr(:)
    integer, optional      , intent(in)    :: logunit
    logical, optional      , intent(in)    :: masterproc

    ! local variables
    integer                         :: rc
    type(dfield_type), pointer      :: dfield_new
    integer                         :: ns, nf
    integer                         :: lsize, num
    integer                         :: status
    integer                         :: fieldCount
    character(cl)                   :: msgstr
    type(ESMF_Field)                :: lfield
    type(ESMF_FieldBundle)          :: fldbun
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=*), parameter :: subname='(dfield_add_strmfld)'
    ! ----------------------------------------------

    allocate(dfield_new, stat=status)
    write(msgstr,*)'allocation error ',__LINE__,':',__FILE__
    if (status /= 0) call shr_sys_abort(msgstr)
    dfield_new%next => dfields
    dfields => dfield_new

    ! determine local size
    lsize = sdat%lsize

    ! initialize dfield_new values
    dfield_new%sdat_stream_index = iunset  ! stream index
    dfield_new%sdat_fldbun_index = iunset  ! index into field bundle for stream(stream_index)
    dfield_new%stream_data1d => null()     ! pointer to stream data 

    ! loop over all input streams and ! determine if the strm_fld is in the attribute vector of stream ns
    do ns = 1, sdat%nstreams
       ! determine if the strm_fld is in the fb_model of stream ns
       call ESMF_FieldBundleGet(sdat%fldbun_model(ns), fieldCount=fieldCount, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) &
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
       allocate(lfieldnamelist(fieldCount))
       call ESMF_FieldBundleGet(sdat%fldbun_model(ns), fieldNameList=lfieldnamelist, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) &
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do nf = 1,fieldcount
          if (trim(strm_fld) == trim(lfieldnamelist(nf))) then
             dfield_new%sdat_stream_index = ns
             dfield_new%sdat_fldbun_index = nf
             allocate(dfield_new%stream_data1d(lsize), stat=status)
             write(msgstr,*)'allocation error ',__LINE__,':',__FILE__
             if (status /= 0) call shr_sys_abort(msgstr)
             strm_ptr => dfield_new%stream_data1d
             if (present(logunit) .and. present(masterproc)) then
                if (masterproc) then
                   write(logunit,*)'(dshr_addfield_add) allocating memory for stream field strm_'//trim(strm_fld)
                end if
             end if
             exit
          end if
       end do
       deallocate(lfieldnamelist)
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
    type(dfield_type), pointer      :: dfield_new
    integer                         :: ns, nf, lsize
    integer                         :: status
    character(cl)                   :: msgstr
    integer                         :: fieldcount
    type(ESMF_Field)                :: lfield
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=*), parameter :: subname='(dfield_add_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    allocate(dfield_new, stat=status)
    write(msgstr,*)'allocation error ',__LINE__,':',__FILE__
    if (status /= 0) call shr_sys_abort(msgstr)
    dfield_new%next => dfields
    dfields => dfield_new

    ! determine local size
    lsize = sdat%lsize

    ! Initialize strm_ptr and state_ptr if it is present.
    ! These will be set to valid values if the relevant fields are
    ! found strm_ptr points to the data in the stream field bundle
    ! that has been spatially and time interpolated to the model mesh

    if (present(strm_ptr )) strm_ptr => null()
    if (present(state_ptr)) state_ptr => null()
    dfield_new%sdat_stream_index = iunset
    dfield_new%sdat_fldbun_index = iunset

    ! always allocate memory for stream_data and initialize it to 0
    ! note that if the attribute vector is not found in the streams

    ! loop over all input streams and ! determine if the strm_fld is in the attribute vector of stream ns
    do ns = 1, sdat%nstreams
       call ESMF_FieldBundleGet(sdat%fldbun_model(ns), fieldCount=fieldCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(lfieldnamelist(fieldCount))
       call ESMF_FieldBundleGet(sdat%fldbun_model(ns), fieldNameList=lfieldnamelist, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do nf = 1,fieldcount
          if (trim(strm_fld) == trim(lfieldnamelist(nf))) then

             ! if strm_fld is in the field bundle of stream ns then
             ! set the field index of the field with the name strm_fld
             ! and set the index of the stream
             dfield_new%sdat_fldbun_index = nf
             dfield_new%sdat_stream_index = ns

             ! set the pointer to the stream data
             allocate(dfield_new%stream_data1d(lsize), stat=status)
             write(msgstr,*)'allocation error ',__LINE__,':',__FILE__
             if (status /= 0) call shr_sys_abort(msgstr)
             dfield_new%stream_data1d(:) = 0._r8
             if (present(strm_ptr)) then
                strm_ptr => dfield_new%stream_data1d
             end if

             ! write output
             if (present(logunit) .and. present(masterproc)) then
                if (masterproc) then
                   write(logunit,*)'(dshr_addfield_add) allocating memory for stream field strm_'//trim(strm_fld)
                end if
             end if
             exit

          end if
       end do
       deallocate(lfieldnamelist)
    end do

    ! Set export state array pointer
    call dshr_state_getfldptr(State, fldname=trim(state_fld), fldptr1=dfield_new%state_data1d, rc=rc)
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
    type(dfield_type), pointer      :: dfield_new
    integer                         :: n, i, ns, nf
    integer                         :: nflds, lsize, num
    integer                         :: status
    character(cl)                   :: msgstr
    integer                         :: fieldcount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    logical                         :: ispresent
    character(len=*), parameter :: subname='(dfield_add_2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    allocate(dfield_new, stat=status)
    write(msgstr,*)'allocation error ',__LINE__,':',__FILE__
    if (status /= 0) call shr_sys_abort(msgstr)
    dfield_new%next => dfields
    dfields => dfield_new

    ! determine local size
    lsize = sdat%lsize

    ! determine stream fldnames array
    nflds = size(strm_flds)

    dfield_new%stream_nflds = nflds
    allocate(dfield_new%stream_fldnames(nflds), stat=status)
    write(msgstr,*)'allocation error ',__LINE__,':',__FILE__
    do n = 1,nflds
       dfield_new%stream_fldnames(n) = strm_flds(n)
    end do

    allocate(dfield_new%sdat_stream_indices(nflds), stat=status)
    write(msgstr,*)'allocation error ',__LINE__,':',__FILE__
    if (status /= 0) call shr_sys_abort(msgstr)

    allocate(dfield_new%sdat_fldbun_indices(nflds), stat=status)
    write(msgstr,*)'allocation error ',__LINE__,':',__FILE__
    if (status /= 0) call shr_sys_abort(msgstr)

    ! loop over all input streams and determine if the strm_flds name
    ! is in the field bundle of stream ns

    ! loop through the field names in strm_flds
    do nf = 1, nflds

       ! loop through input streams
       do ns = 1, sdat%nstreams

          ! determine which stream the field with name dfield%stream_fldnames(nf) is in
          call ESMF_FieldBundleGet(sdat%fldbun_model(ns), fieldName=trim(dfield_new%stream_fldnames(nf)), &
               isPresent=isPresent, rc=rc)
          if (ispresent) then
             ! if field is present in stream - determine the index in the field bundle of this field
             dfield_new%sdat_stream_indices(nf) = ns
             call ESMF_FieldBundleGet(sdat%fldbun_model(ns), fieldCount=fieldCount, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             allocate(lfieldnamelist(fieldCount))
             call ESMF_FieldBundleGet(sdat%fldbun_model(ns), fieldNameList=lfieldnamelist, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do n = 1,fieldcount
                if (trim(strm_flds(nf)) == trim(lfieldnamelist(n))) then
                   dfield_new%sdat_fldbun_indices(nf) = n
                end if
             end do
             deallocate(lfieldnamelist)
             exit ! go to the next fld
          end if
          if (present(logunit) .and. present(masterproc)) then
             if (masterproc) then
                write(logunit,*)'(dshr_addfield_add) using stream field strm_'//&
                     trim(strm_flds(nf))//' for 2d '//trim(state_fld) 
             end if
          end if
       end do
    end do

    ! Set export state array pointer
    call dshr_state_getfldptr(State, fldname=trim(state_fld), fldptr2=dfield_new%state_data2d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dfield_new%state_data2d(nflds,lsize) = 0._r8
    if (present(logunit) .and. present(masterproc)) then
       if (masterproc) then
          write(logunit,*)'(dshr_addfield_add) setting pointer for export state '//trim(state_fld)
       end if
    end if

    ! Return array pointer if argument is present
    if (present(state_ptr)) then
       state_ptr => dfield_new%state_data2d
    end if

  end subroutine dshr_dfield_add_2d

  !===============================================================================

  subroutine dshr_dfield_copy(dfields, sdat, rc)

    ! Copy stream data into dfield data type for each element of dfields
    ! This routine will do one of the following
    ! - populate the export state data (dfield%state_data1d or dfield%state_data2d) 
    !   with the stream field data
    ! - populate the dfield stream field (dfield%stream_data1d) with spatially and
    !   time interpolate stream field data

    ! input/output variables
    type(dfield_type)      , pointer     :: dfields
    type(shr_strdata_type) , intent(in)  :: sdat
    integer                , intent(out) :: rc

    ! local variables
    type(ESMF_field)           :: lfield
    type(dfield_type), pointer :: dfield
    real(r8), pointer          :: data1d(:)
    integer                    :: n, nf, i, k
    integer                    :: fldbun_index
    integer                    :: stream_index
    logical                    :: ispresent
    integer                    :: ns
    character(len=CL)          :: msgstr          ! temporary
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Loop over all dfield entries and fill in stream_data and state_data1d or state_data2d arrays
    dfield => dfields ! note that dfields is the head of the linked list
    do while (associated(dfield))

       ! Map the stream data to the state data
       if (associated(dfield%state_data1d)) then
          stream_index = dfield%sdat_stream_index
          fldbun_index = dfield%sdat_fldbun_index
          if (stream_index /= iunset .and. fldbun_index /= iunset) then
             call dshr_fldbun_getfieldn(sdat%fldbun_model(stream_index), fldbun_index, lfield, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call dshr_field_getfldptr(lfield, fldptr1=data1d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dfield%stream_data1d(:) = data1d(:)
             dfield%state_data1d(:) = data1d(:)
          end if
       else if (associated(dfield%state_data2d)) then
          do nf = 1,dfield%stream_nflds
             stream_index = dfield%sdat_stream_indices(nf)
             fldbun_index = dfield%sdat_fldbun_indices(nf)
             call dshr_fldbun_getfieldn(sdat%fldbun_model(stream_index), fldbun_index, lfield, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call dshr_field_getfldptr(lfield, fldptr1=data1d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dfield%state_data2d(nf,:) = data1d(:)
          end do

       else if (associated(dfield%stream_data1d)) then
          stream_index = dfield%sdat_stream_index
          fldbun_index = dfield%sdat_fldbun_index
          if (stream_index /= iunset .and. fldbun_index /= iunset) then
             call dshr_fldbun_getfieldn(sdat%fldbun_model(stream_index), fldbun_index, lfield, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call dshr_field_getfldptr(lfield, fldptr1=data1d, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dfield%stream_data1d(:) = data1d(:)
          end if
       end if
       dfield => dfield%next

    end do

  end subroutine dshr_dfield_copy

end module dshr_dfield_mod
