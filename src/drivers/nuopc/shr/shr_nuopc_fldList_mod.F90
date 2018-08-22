module shr_nuopc_fldList_mod
  use shr_kind_mod   , only : CX => shr_kind_CX, CS=>shr_kind_CS, CL=>shr_kind_cl

  implicit none
  private

  integer, parameter, public :: CSS = 256  ! use longer short character
  integer, parameter, public :: CLL = 1024

  public :: shr_nuopc_fldList_Concat
  public :: shr_nuopc_fldList_Realize
  public :: shr_nuopc_fldList_AddFld
  public :: shr_nuopc_fldList_AddDomain
  public :: shr_nuopc_fldList_AddMetadata
  public :: shr_nuopc_fldList_GetMetadata
  public :: shr_nuopc_fldList_AddMap
  public :: shr_nuopc_fldList_Deactivate
  public :: shr_nuopc_fldList_GetFldNames
  public :: shr_nuopc_fldList_GetNumFlds
  public :: shr_nuopc_fldList_GetFldInfo

  !-----------------------------------------------
  ! Metadata array
  !-----------------------------------------------

  character(len=*), parameter  :: undef     = 'undefined'
  integer         , parameter  :: nmax      = 1000        ! maximum number of entries in metadta_entry
  integer                      :: n_entries = 0           ! actual number of entries in metadta_entry
  character(len=CSS)           :: shr_nuopc_fldList_Metadata(nmax,4) = undef

  !-----------------------------------------------
  ! Maximum number of components, mappers
  !-----------------------------------------------

  integer          , public, parameter :: ncomps_max = 8
  integer          , public, parameter :: mapunset=0
  integer          , public, parameter :: nmappers=6
  integer          , public, parameter :: mapbilnr=1
  integer          , public, parameter :: mapconsf=2
  integer          , public, parameter :: mapconsd=3
  integer          , public, parameter :: mappatch=4
  integer          , public, parameter :: mapfcopy=5
  integer          , public, parameter :: mapfiler=6
  character(len=*) , public, parameter :: mapnames(nmappers) = (/'bilnr','consf','consd','patch','fcopy','filer'/)

  !-----------------------------------------------
  ! Types and instantiations that determine fields, mappings, mergings
  !-----------------------------------------------

  type shr_nuopc_fldList_entry_type
     character(CS) :: stdname
     character(CS) :: shortname
     logical       :: active = .true.
     ! Mapping fldsFr data - for mediator import fields
     integer       :: mapindex(ncomps_max) = mapunset
     character(CS) :: mapnorm(ncomps_max) = 'unset'
     character(CX) :: mapfile(ncomps_max) = 'unset'
     ! Merging fldsTo data - for mediator export fields
     character(CX) :: merge_fields(ncomps_max)    = 'unset'
     character(CS) :: merge_types(ncomps_max)     = 'unset'
     character(CS) :: merge_fracnames(ncomps_max) = 'unset'
  end type shr_nuopc_fldList_entry_type
  public :: shr_nuopc_fldList_entry_type

  ! The above would be the field name to merge from
  ! e.g. for Sa_z in lnd
  !    merge_field(compatm) = 'Sa_z'
  !    merge_type(comptm) = 'copy'  (could also have 'copy_with_weighting')

  type shr_nuopc_fldList_type
     type (shr_nuopc_fldList_entry_type), pointer :: flds(:)
  end type shr_nuopc_fldList_type
  public :: shr_nuopc_fldList_type

  interface shr_nuopc_fldList_GetFldInfo ; module procedure &
       shr_nuopc_fldList_GetFldInfo_general, &
       shr_nuopc_fldList_GetFldInfo_stdname, &
       shr_nuopc_fldList_GetFldInfo_merging
  end interface

  integer           :: dbrc
  character(len=CL) :: infostr
  character(len=*),parameter :: u_FILE_u = &
     __FILE__

!================================================================================
contains
!================================================================================

  subroutine shr_nuopc_fldList_Concat(fldsFr, fldsTo, concat_src, concat_dst, flds_scalar_name)
    ! Returns new concatentated colon delimited field lists
    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_ERROR
    ! input/output parameters:
    type(shr_nuopc_fldList_type) , intent(in)    :: fldsFr
    type(shr_nuopc_fldList_type) , intent(in)    :: fldsTo
    character(len=*)             , intent(in)    :: flds_scalar_name
    character(len=*)             , intent(inout) :: concat_src
    character(len=*)             , intent(inout) :: concat_dst

    ! local variables
    integer :: n
    character(len=*),parameter :: subname = '(shr_nuopc_fldList_concat) '
    !-------------------------------------------------------------------------------

    do n = 1,size(FldsFr%flds)
       if (trim(FldsFr%flds(n)%shortname) /= flds_scalar_name) then
          if (len_trim(concat_src) + len_trim(FldsFr%flds(n)%shortname) + 1 >= len(concat_src)) then
             call ESMF_LogWrite(subname//': ERROR: max len of fldlist has been exceeded', ESMF_LOGMSG_ERROR, line=__LINE__, file= u_FILE_u, rc=dbrc)
             return
          end if
          if (trim(concat_src) == '') then
             concat_src = trim(FldsFr%flds(n)%shortname)
          else
             concat_src = trim(concat_src)//':'//trim(FldsFr%flds(n)%shortname)
          end if
       end if
    end do

    do n = 1,size(FldsTo%flds)
       if (trim(FldsTo%flds(n)%shortname) /= flds_scalar_name) then
          if (len_trim(concat_dst) + len_trim(FldsTo%flds(n)%shortname) + 1 >= len(concat_dst)) then
             call ESMF_LogWrite(subname//': ERROR: max len of fldlist has been exceeded', &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file= u_FILE_u, rc=dbrc)
             return
          end if
          if (trim(concat_dst) == '') then
             concat_dst = trim(FldsTo%flds(n)%shortname)
          else
             concat_dst = trim(concat_dst)//':'//trim(FldsTo%flds(n)%shortname)
          end if
       end if
    end do

  end subroutine shr_nuopc_fldList_Concat

  !===============================================================================

  subroutine shr_nuopc_fldList_AddDomain(fldlist, fldname, longname, stdname, units)

    ! Returns new concatentated field and map lists
    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_ERROR
    ! input/output parameters:
    character(len=*),intent(inout)       :: fldlist   ! output field name
    character(len=*),intent(in)          :: fldname   ! fldname to add to fldlist
    character(len=*),intent(in),optional :: longname
    character(len=*),intent(in),optional :: stdname
    character(len=*),intent(in),optional :: units

    ! local variables
    character(len=*),parameter :: subname = '(fldList_AddDomain) '
    !-------------------------------------------------------------------------------

    if (len_trim(fldlist) + len_trim(fldname) + 1 >= len(fldlist)) then
       call ESMF_LogWrite(subname//': ERROR: max len of fldlist has been exceeded', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file= u_FILE_u, rc=dbrc)
       return
    end if

    if (trim(fldlist) == '') then
       fldlist = trim(fldname)
    else
       fldlist = trim(fldlist)//':'//trim(fldname)
    end if

    if (present(longname) .and. present(stdname) .and. present(units)) then
       call shr_nuopc_fldList_AddMetadata(trim(fldname), longname, stdname, units)
    endif

  end subroutine shr_nuopc_fldList_AddDomain

  !===============================================================================

  subroutine shr_nuopc_fldList_AddMetadata(fldname , longname, stdname, units)

    use NUOPC, only : NUOPC_FieldDictionaryAddEntry, NUOPC_FieldDictionaryHasEntry
    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
    use ESMF, only : ESMF_LOGMSG_ERROR, ESMF_FAILURE
    ! input/output parameters:
    character(len=*), intent(in) :: fldname
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: stdname
    character(len=*), intent(in) :: units

    ! local variables
    integer :: n
    logical :: found,FDfound
    integer :: rc
    character(len=*),parameter :: subname = '(fldList_AddMetadata) '
    !-------------------------------------------------------------------------------

    FDfound = .true.
    if (.not.NUOPC_FieldDictionaryHasEntry(fldname)) then
       FDfound = .false.
       call ESMF_LogWrite(subname//': Add:'//trim(fldname), ESMF_LOGMSG_INFO, rc=dbrc)
       call NUOPC_FieldDictionaryAddEntry(standardName=fldname, canonicalUnits=units, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    endif

    found = .false.
    ! only do the search if it was already in the FD.  If it wasn't,
    ! then assume it's also not in the metadata table.
    if (FDfound) then
       n = 1
       do while (n <= n_entries .and. .not.found)
          if (fldname == shr_nuopc_fldList_Metadata(n,1)) found=.true.
          n = n + 1
       enddo
    endif

    if (.not. found) then
       n_entries = n_entries + 1
       if (n_entries > nmax) then
          write(infostr,*) subname,' ERROR: n_entries= ',n_entries,' nmax = ',nmax,' fldname= ',trim(fldname)
          call ESMF_LogWrite(trim(infostr),ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          write(infostr,*) subname,' ERROR: n_entries gt nmax'
          call ESMF_LogWrite(trim(infostr),ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
       end if
       shr_nuopc_fldList_Metadata(n_entries,1) = trim(fldname)
       shr_nuopc_fldList_Metadata(n_entries,2) = trim(longname)
       shr_nuopc_fldList_Metadata(n_entries,3) = trim(stdname )
       shr_nuopc_fldList_Metadata(n_entries,4) = trim(units   )
    endif

  end subroutine shr_nuopc_fldList_AddMetadata

  !===============================================================================

  subroutine shr_nuopc_fldList_GetMetadata(shortname, longname, stdname, units)

    ! !USES:
    use shr_string_mod , only : shr_string_lastindex
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*), intent(in)  :: shortname
    character(len=*),optional, intent(out) :: longname
    character(len=*),optional, intent(out) :: stdname
    character(len=*),optional, intent(out) :: units

    !EOP

    !--- local ---
    integer :: i,n
    character(len=CSS) :: llongname, lstdname, lunits, lshortname  ! local copies
    character(len=*),parameter :: unknown = 'unknown'
    logical :: found
    character(len=*),parameter :: subname = '(shr_nuopc_fldList_GetMetadata) '

    !--- define field metadata (name, long_name, standard_name, units) ---

    llongname = trim(unknown)
    lstdname  = trim(unknown)
    lunits    = trim(unknown)

    found = .false.

    if (.not.found) then
       i = 1
       do while (i <= n_entries .and. .not.found)
          lshortname = trim(shortname)
          if (trim(lshortname) == trim(shr_nuopc_fldList_Metadata(i,1))) then
             llongname = trim(shr_nuopc_fldList_Metadata(i,2))
             lstdname  = trim(shr_nuopc_fldList_Metadata(i,3))
             lunits    = trim(shr_nuopc_fldList_Metadata(i,4))
             found     =.true.
          end if
          i = i + 1
       end do
    endif

    if (.not.found) then
       i = 1
       do while (i <= n_entries .and. .not.found)
          n = shr_string_lastIndex(shortname,"_")
          lshortname = ""
          if (n < len_trim(shortname)) lshortname = shortname(n+1:len_trim(shortname))
          if (trim(lshortname) == trim(shr_nuopc_fldList_Metadata(i,1))) then
             llongname = trim(shr_nuopc_fldList_Metadata(i,2))
             lstdname  = trim(shr_nuopc_fldList_Metadata(i,3))
             lunits    = trim(shr_nuopc_fldList_Metadata(i,4))
             found     = .true.
          end if
          i = i + 1
       end do
    endif

    if (present(longname)) then
       longname = trim(llongname)
    endif
    if (present(stdname))  then
       stdname = trim(lstdname)
    endif
    if (present(units)) then
       units = trim(lunits)
    endif

  end subroutine shr_nuopc_fldList_GetMetadata

  !================================================================================

  subroutine shr_nuopc_fldList_AddFld(flds, stdname, shortname, fldindex, &
       merge_from1, merge_field1, merge_type1, merge_fracname1, &
       merge_from2, merge_field2, merge_type2, merge_fracname2, &
       merge_from3, merge_field3, merge_type3, merge_fracname3, &
       merge_from4, merge_field4, merge_type4, merge_fracname4)

    ! ----------------------------------------------
    ! Add an entry to to the flds array
    ! Use pointers to create an extensible allocatable array.
    ! to allow the size of flds to grow, the process for
    ! adding a new field is:
    ! 1) allocate newflds to be N (one element larger than flds)
    ! 2) copy flds into first N-1 elements of newflds
    ! 3) newest flds entry is Nth element of newflds
    ! 4) deallocate / nullify flds
    ! 5) point flds => newflds
    ! ----------------------------------------------

    type(shr_nuopc_fldList_entry_type) , pointer                :: flds(:)
    character(len=*)                   , intent(in)             :: stdname
    character(len=*)                   , intent(in)  , optional :: shortname
    integer                            , intent(out) , optional :: fldindex
    integer                            , intent(in)  , optional :: merge_from1
    character(len=*)                   , intent(in)  , optional :: merge_field1
    character(len=*)                   , intent(in)  , optional :: merge_type1
    character(len=*)                   , intent(in)  , optional :: merge_fracname1
    integer                            , intent(in)  , optional :: merge_from2
    character(len=*)                   , intent(in)  , optional :: merge_field2
    character(len=*)                   , intent(in)  , optional :: merge_type2
    character(len=*)                   , intent(in)  , optional :: merge_fracname2
    integer                            , intent(in)  , optional :: merge_from3
    character(len=*)                   , intent(in)  , optional :: merge_field3
    character(len=*)                   , intent(in)  , optional :: merge_type3
    character(len=*)                   , intent(in)  , optional :: merge_fracname3
    integer                            , intent(in)  , optional :: merge_from4
    character(len=*)                   , intent(in)  , optional :: merge_field4
    character(len=*)                   , intent(in)  , optional :: merge_type4
    character(len=*)                   , intent(in)  , optional :: merge_fracname4

    ! local variables
    integer :: n,oldsize,id
    type(shr_nuopc_fldList_entry_type), pointer :: newflds(:)
    character(len=*), parameter :: subname='(fldList_AddFld)'
    ! ----------------------------------------------

    if (associated(flds)) then
       oldsize = size(flds)
    else
       oldsize = 0
    end if
    id = oldsize + 1

    ! 1) allocate newfld to be size (one element larger than input flds)
    allocate(newflds(id))

    ! 2) copy flds into first N-1 elements of newflds
    do n = 1,oldsize
       newflds(n)%stdname            = flds(n)%stdname
       newflds(n)%shortname          = flds(n)%shortname
       newflds(n)%active             = flds(n)%active
       newflds(n)%mapindex(:)        = flds(n)%mapindex(:)
       newflds(n)%mapnorm(:)         = flds(n)%mapnorm(:)
       newflds(n)%mapfile(:)         = flds(n)%mapfile(:)
       newflds(n)%merge_fields(:)    = flds(n)%merge_fields(:)
       newflds(n)%merge_types(:)     = flds(n)%merge_types(:)
       newflds(n)%merge_fracnames(:) = flds(n)%merge_fracnames(:)
    end do

    ! 3) deallocate / nullify flds
    if (oldsize >  0) then
       deallocate(flds)
       nullify(flds)
    end if

    ! 4) point flds => new_flds
    flds => newflds

    ! 5) now update flds information for new entry
    flds(id)%stdname   = trim(stdname)
    if (present(shortname)) then
       flds(id)%shortname = trim(shortname)
    else
       flds(id)%shortname = trim(stdname)
    end if
    if (present(fldindex)) then
       fldindex = id
    end if
    if (present(merge_from1) .and. present(merge_field1) .and. present(merge_type1)) then
       n = merge_from1
       flds(id)%merge_fields(n) = merge_field1
       flds(id)%merge_types(n) = merge_type1
       if (present(merge_fracname1)) then
          flds(id)%merge_fracnames(n) = merge_fracname1
       end if
    end if
    if (present(merge_from2) .and. present(merge_field2) .and. present(merge_type2)) then
       n = merge_from2
       flds(id)%merge_fields(n) = merge_field2
       flds(id)%merge_types(n) = merge_type2
       if (present(merge_fracname2)) then
          flds(id)%merge_fracnames(n) = merge_fracname2
       end if
    end if
    if (present(merge_from3) .and. present(merge_field3) .and. present(merge_type3)) then
       n = merge_from3
       flds(id)%merge_fields(n) = merge_field3
       flds(id)%merge_types(n) = merge_type3
       if (present(merge_fracname3)) then
          flds(id)%merge_fracnames(n) = merge_fracname3
       end if
    end if
    if (present(merge_from4) .and. present(merge_field4) .and. present(merge_type4)) then
       n = merge_from4
       flds(id)%merge_fields(n) = merge_field4
       flds(id)%merge_types(n) = merge_type4
       if (present(merge_fracname4)) then
          flds(id)%merge_fracnames(n) = merge_fracname4
       end if
    end if
  end subroutine shr_nuopc_fldList_AddFld

  !================================================================================

  subroutine shr_nuopc_fldList_AddMap(fld, srccomp, destcomp, mapindex, mapnorm, mapfile)
    type(shr_nuopc_fldList_entry_type) , intent(inout) :: fld
    integer                            , intent(in)    :: srccomp
    integer                            , intent(in)    :: destcomp
    integer                            , intent(in)    :: mapindex
    character(len=*)                   , intent(in)    :: mapnorm
    character(len=*)                   , intent(in)    :: mapfile

    ! local variables
    logical :: mapset
    character(len=*),parameter  :: subname='(fldList_AddMap)'
    ! ----------------------------------------------

    ! Note - default values are already set for the fld entries - so only non-default
    ! values need to be set below
    ! If mapindex is mapfcopy - create a redistribution route handle
    ! If mapfile is idmap - create a redistribution route nhandle
    ! If mapfile is unset then create the mapping route handle at run time

    fld%mapindex(destcomp) = mapindex
    fld%mapfile(destcomp)  = trim(mapfile)
    fld%mapnorm(destcomp)  = trim(mapnorm)

    ! overwrite values if appropriate
    if (fld%mapindex(destcomp) == mapfcopy) then
       fld%mapfile(destcomp) = 'unset'
       fld%mapnorm(destcomp) = 'unset'
    else if (trim(fld%mapfile(destcomp)) == 'idmap') then
       fld%mapindex(destcomp) = mapfcopy
       fld%mapnorm(destcomp) = 'unset'
    end if
  end subroutine shr_nuopc_fldList_AddMap

  !================================================================================

  subroutine shr_nuopc_fldList_Realize(state, fldList, flds_scalar_name, flds_scalar_num, &
       grid, mesh, tag, rc)
    use NUOPC, only : NUOPC_GetStateMemberLists, NUOPC_IsConnected, NUOPC_Realize
    use NUOPC, only : NUOPC_GetAttribute
    use ESMF, only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF, only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Grid, ESMF_Mesh
    use ESMF, only : ESMF_StateGet, ESMF_LogFoundError
    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_FAILURE, ESMF_LOGERR_PASSTHRU
    use ESMF, only : ESMF_LOGMSG_INFO, ESMF_StateRemove, ESMF_SUCCESS
    use med_constants_mod       , only : dbug_flag=>med_constants_dbug_flag
    
    type(ESMF_State)            , intent(inout)            :: state
    type(shr_nuopc_fldlist_type), intent(in)               :: fldList
    character(len=*)            , intent(in)               :: flds_scalar_name
    integer                     , intent(in)               :: flds_scalar_num
    character(len=*)            , intent(in)               :: tag
    integer                     , intent(inout)            :: rc
    type(ESMF_Grid)             , intent(in)    , optional :: grid
    type(ESMF_Mesh)             , intent(in)    , optional :: mesh

    ! local variables
    integer                         :: n, nflds
    integer                         :: itemCount
    type(ESMF_Field)                :: field
    character(CS)                   :: shortname
    character(CS)                   :: stdname
    character(ESMF_MAXSTR)          :: transferAction
    character(ESMF_MAXSTR), pointer :: StandardNameList(:)
    character(ESMF_MAXSTR), pointer :: ConnectedList(:)
    character(ESMF_MAXSTR), pointer :: NameSpaceList(:)
    character(ESMF_MAXSTR), pointer :: itemNameList(:)
    character(len=*),parameter  :: subname='(shr_nuopc_fldList_Realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (present(grid) .and. present(mesh)) then
       call ESMF_LogWrite(trim(subname)//trim(tag)//": ERROR both grid and mesh not allowed", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
       rc = ESMF_FAILURE
       return
    endif

    nullify(StandardNameList)
    nullify(ConnectedList)
    nullify(NameSpaceList)
    nullify(ItemNameList)

    call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    write(infostr,'(i6)') itemCount
    call ESMF_LogWrite(trim(subname)//trim(tag)//" count = "//trim(infostr), ESMF_LOGMSG_INFO, rc=dbrc)
    if (itemCount > 0) then
       allocate(itemNameList(itemCount))
       call ESMF_StateGet(state, itemNameList=itemNameList, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       do n = 1,itemCount
          call ESMF_LogWrite(trim(subname)//trim(tag)//" itemNameList = "//trim(itemNameList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
       enddo
       deallocate(itemNameList)
    endif

#if (1 == 0)
    call NUOPC_GetStateMemberLists(state, StandardNameList=StandardNameList, ConnectedList=ConnectedList, &
         NamespaceList=NamespaceList, itemNameList=itemNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    write(infostr,'(i6)') size(StandardNameList)
    call ESMF_LogWrite(trim(subname)//trim(tag)//" size = "//trim(infostr), ESMF_LOGMSG_INFO, rc=dbrc)

    do n = 1,size(StandardNameList)
       call ESMF_LogWrite(trim(subname)//trim(tag)//" StandardNameList = "//trim(StandardNameList(n)), &
            ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    do n = 1,size(ConnectedList)
       call ESMF_LogWrite(trim(subname)//trim(tag)//" ConnectedList = "//trim(ConnectedList(n)), &
            ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    do n = 1,size(NamespaceList)
       call ESMF_LogWrite(trim(subname)//trim(tag)//" NamespaceList = "//trim(NamespaceList(n)), &
            ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    do n = 1,size(ItemnameList)
       call ESMF_LogWrite(trim(subname)//trim(tag)//" ItemnameList = "//trim(ItemnameList(n)), &
            ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
#endif

    nflds = size(fldList%flds)

    do n = 1, nflds
       shortname = fldList%flds(n)%shortname

       if (fldList%flds(n)%active) then
          ! call ESMF_LogWrite(subname//' fld = '//trim(shortname), ESMF_LOGMSG_INFO, rc=dbrc)

          if (NUOPC_IsConnected(state, fieldName=shortname)) then

             call ESMF_StateGet(state, field=field, itemName=trim(shortname), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

             call NUOPC_GetAttribute(field, name="TransferActionGeomObject", value=transferAction, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

             if (trim(transferAction) == "accept") then  ! accept

                call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(shortname)//" is connected, grid/mesh TBD", &
                     ESMF_LOGMSG_INFO, rc=dbrc)

             else   ! provide

                if (shortname == trim(flds_scalar_name)) then
                   call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(shortname)//" is connected on root pe", &
                        ESMF_LOGMSG_INFO, rc=dbrc)
                   call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
                   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
                elseif (present(grid)) then
                   call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(shortname)//" is connected using grid", &
                        ESMF_LOGMSG_INFO, rc=dbrc)
                   ! Create the field
                   field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=shortname,rc=rc)
                   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
                elseif (present(mesh)) then
                   call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(shortname)//" is connected using mesh", &
                        ESMF_LOGMSG_INFO, rc=dbrc)
                   ! Create the field
                   field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=shortname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
                else
                   call ESMF_LogWrite(trim(subname)//trim(tag)//": ERROR grid or mesh expected", &
                        ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
                   rc = ESMF_FAILURE
                   return
                endif

                ! NOW call NUOPC_Realize
                call NUOPC_Realize(state, field=field, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

                ! call ESMF_FieldPrint(field=field, rc=rc)
                ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

             endif

          else

             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(shortname) // " is not connected.", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             call ESMF_StateRemove(state, (/shortname/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

          end if

       end if
    end do

    call ESMF_LogWrite(subname//' done ', ESMF_LOGMSG_INFO, rc=dbrc)

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8
      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), &
           grid=grid, &
           typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), &
           ungriddedUBound=(/flds_scalar_num/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine shr_nuopc_fldList_Realize

  !================================================================================

  subroutine shr_nuopc_fldList_Deactivate(fldList, flds_scalar_name)
    ! ----------------------------------------------
    ! set active flag to .false. for all fields other than flds_scalar_name
    ! ----------------------------------------------
    type(shr_nuopc_fldList_type) , intent(in) :: fldList
    character(len=*)             , intent(in) :: flds_scalar_name

    ! local variables
    integer :: n
    character(len=*), parameter :: subname='(shr_nuopc_fldList_Deactivate)'
    ! ----------------------------------------------

    do n = 1,size(fldList%flds)
       if (trim(fldList%flds(n)%shortname) /= flds_scalar_name) then
          fldList%flds(n)%active = .false.
       end if
    end do
  end subroutine shr_nuopc_fldList_Deactivate

  !================================================================================

  subroutine shr_nuopc_fldList_GetFldInfo_general(fldList, fldindex, active, stdname, shortname)
    ! ----------------------------------------------
    ! Get field info
    ! ----------------------------------------------
    type(shr_nuopc_fldList_type) , intent(in)  :: fldList
    integer                      , intent(in)  :: fldindex
    logical                      , intent(out) :: active
    character(len=*)             , intent(out) :: stdname
    character(len=*)             , intent(out) :: shortname

    ! local variables
    character(len=*), parameter :: subname='(shr_nuopc_fldList_GetFldInfo_general)'
    ! ----------------------------------------------

    active    = fldList%flds(fldindex)%active
    stdname   = fldList%flds(fldindex)%stdname
    shortname = fldList%flds(fldindex)%shortname
  end subroutine shr_nuopc_fldList_GetFldInfo_general

  subroutine shr_nuopc_fldList_GetFldInfo_stdname(fldList, fldindex, stdname)
    ! ----------------------------------------------
    ! Get field info
    ! ----------------------------------------------
    type(shr_nuopc_fldList_type) , intent(in)  :: fldList
    integer                      , intent(in)  :: fldindex
    character(len=*)             , intent(out) :: stdname

    ! local variables
    character(len=*), parameter :: subname='(shr_nuopc_fldList_GetFldInfo_stdname)'
    ! ----------------------------------------------

    stdname   = fldList%flds(fldindex)%stdname
  end subroutine shr_nuopc_fldList_GetFldInfo_stdname

  subroutine shr_nuopc_fldList_GetFldInfo_merging(fldList, fldindex, compsrc, merge_field, merge_type, merge_fracname)
    ! ----------------------------------------------
    ! Get field merge info
    ! ----------------------------------------------
    type(shr_nuopc_fldList_type) , intent(in)  :: fldList
    integer                      , intent(in)  :: fldindex
    integer                      , intent(in)  :: compsrc
    character(len=*)             , intent(out) :: merge_field
    character(len=*)             , intent(out) :: merge_type
    character(len=*)             , intent(out) :: merge_fracname

    ! local variables
    character(len=*), parameter :: subname='(shr_nuopc_fldList_GetFldInfo_merging)'
    ! ----------------------------------------------

    merge_field    = fldList%flds(fldindex)%merge_fields(compsrc)
    merge_type     = fldList%flds(fldindex)%merge_types(compsrc)
    merge_fracname = fldList%flds(fldindex)%merge_fracnames(compsrc)
  end subroutine shr_nuopc_fldList_GetFldInfo_merging

  !================================================================================

  integer function shr_nuopc_fldList_GetNumFlds(fldList)
    ! ----------------------------------------------
    ! Get number of fields
    ! ----------------------------------------------
    type(shr_nuopc_fldList_type), intent(in)  :: fldList
    shr_nuopc_fldList_GetNumFlds = size(fldList%flds)
  end function shr_nuopc_fldList_GetNumFlds

  !================================================================================

  subroutine shr_nuopc_fldList_GetFldNames(flds, fldnames)
    type(shr_nuopc_fldList_entry_type), intent(in) :: flds(:)
    character(len=*), pointer :: fldnames(:)
    integer :: n
    do n = 1,size(flds)
       fldnames(n) = trim(flds(n)%shortname)
    end do
  end subroutine shr_nuopc_fldList_GetFldNames


end module shr_nuopc_fldList_mod
