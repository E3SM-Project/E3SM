module esmflds

  use shr_kind_mod , only : CX => shr_kind_CX, CS=>shr_kind_CS, CL=>shr_kind_cl
  use shr_sys_mod  , only : shr_sys_abort

  implicit none
  private

  !-----------------------------------------------
  ! Set components
  !-----------------------------------------------

  integer, public, parameter  :: ncomps  = 8
  integer, public, parameter  :: compmed = 1
  integer, public, parameter  :: compatm = 2
  integer, public, parameter  :: complnd = 3
  integer, public, parameter  :: compocn = 4
  integer, public, parameter  :: compice = 5
  integer, public, parameter  :: comprof = 6
  integer, public, parameter  :: compwav = 7
  integer, public, parameter  :: compglc = 8

  character(len=*), public, parameter :: compname(ncomps) = &
       (/'med','atm','lnd','ocn','ice','rof','wav','glc'/)

  !-----------------------------------------------
  ! Set mappers
  !-----------------------------------------------

  integer , public, parameter :: mapunset       = 0
  integer , public, parameter :: mapbilnr       = 1
  integer , public, parameter :: mapconsf       = 2
  integer , public, parameter :: mapconsd       = 3
  integer , public, parameter :: mappatch       = 4
  integer , public, parameter :: mapfcopy       = 5
  integer , public, parameter :: mapfiler       = 6
  integer , public, parameter :: mapnstod       = 7 ! nearest source to destination
  integer , public, parameter :: mapnstod_consd = 8 ! nearest source to destination followed by conservative dst
  integer , public, parameter :: mapnstod_consf = 9 ! nearest source to destination followed by conservative frac
  integer , public, parameter :: nmappers       = 9

  character(len=*) , public, parameter :: mapnames(nmappers) = & 
       (/'bilnr', 'consf', 'consd', 'patch', 'fcopy', 'filer', 'nstod', 'nstod_consd', 'nstod_consf'/)

  !-----------------------------------------------
  ! Set coupling mode
  !-----------------------------------------------

  character(len=10), public :: coupling_mode ! valid values are [cesm,nems_orig,nems_frac]

  !-----------------------------------------------
  ! PUblic methods
  !-----------------------------------------------

  public :: shr_nuopc_fldList_AddFld
  public :: shr_nuopc_fldList_AddMap
  public :: shr_nuopc_fldList_AddMrg
  public :: shr_nuopc_fldList_AddMetadata
  public :: shr_nuopc_fldList_GetMetadata
  public :: shr_nuopc_fldList_GetFldNames
  public :: shr_nuopc_fldList_GetNumFlds
  public :: shr_nuopc_fldList_GetFldInfo
  public :: shr_nuopc_fldList_Realize
  public :: shr_nuopc_fldList_Document_Mapping
  public :: shr_nuopc_fldList_Document_Merging

  !-----------------------------------------------
  ! Metadata array
  !-----------------------------------------------

  integer, public , parameter  :: CSS = 256  ! use longer short character
  character(len=*), parameter  :: undef     = 'undefined'
  integer         , parameter  :: nmax      = 1000        ! maximum number of entries in metadta_entry
  integer                      :: n_entries = 0           ! actual number of entries in metadta_entry
  character(len=CSS)           :: shr_nuopc_fldList_Metadata(nmax,4) = undef

  !-----------------------------------------------
  ! Types and instantiations that determine fields, mappings, mergings
  !-----------------------------------------------

  type, public :: shr_nuopc_fldList_entry_type
     character(CS) :: stdname
     character(CS) :: shortname

     ! Mapping fldsFr data - for mediator import fields
     integer       :: mapindex(ncomps) = mapunset
     character(CS) :: mapnorm(ncomps) = 'unset'
     character(CX) :: mapfile(ncomps) = 'unset'

     ! Merging fldsTo data - for mediator export fields
     character(CS) :: merge_fields(ncomps)    = 'unset'
     character(CS) :: merge_types(ncomps)     = 'unset'
     character(CS) :: merge_fracnames(ncomps) = 'unset'
  end type shr_nuopc_fldList_entry_type

  ! The above would be the field name to merge from
  ! e.g. for Sa_z in lnd
  !    merge_field(compatm) = 'Sa_z'
  !    merge_type(comptm) = 'copy'  (could also have 'copy_with_weighting')

  type, public :: shr_nuopc_fldList_type
     type (shr_nuopc_fldList_entry_type), pointer :: flds(:)
  end type shr_nuopc_fldList_type

  interface shr_nuopc_fldList_GetFldInfo ; module procedure &
       shr_nuopc_fldList_GetFldInfo_general, &
       shr_nuopc_fldList_GetFldInfo_stdname, &
       shr_nuopc_fldList_GetFldInfo_merging
  end interface

  !-----------------------------------------------
  ! Instantiate derived types
  !-----------------------------------------------
  type (shr_nuopc_fldList_type), public :: fldListTo(ncomps) ! advertise fields to components
  type (shr_nuopc_fldList_type), public :: fldListFr(ncomps) ! advertise fields from components

  type (shr_nuopc_fldList_type), public :: fldListMed_aoflux 
  type (shr_nuopc_fldList_type), public :: fldListMed_ocnalb 

  integer                    :: dbrc
  character(len=CL)          :: infostr
  character(len=*),parameter :: u_FILE_u = &
     __FILE__

!================================================================================
contains
!================================================================================

  subroutine shr_nuopc_fldList_AddMetadata(fldname , longname, stdname, units)

    use NUOPC , only : NUOPC_FieldDictionaryAddEntry, NUOPC_FieldDictionaryHasEntry
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
    use ESMF  , only : ESMF_LOGMSG_ERROR, ESMF_FAILURE

    ! input/output parameters:
    character(len=*), intent(in) :: fldname
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: stdname
    character(len=*), intent(in) :: units

    ! local variables
    integer :: n
    logical :: found,FDfound
    integer :: rc
    character(len=*),parameter :: subname = '(shr_nuopc_fldList_AddMetadata) '
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

    use shr_string_mod , only : shr_string_lastindex

    ! input/output variables
    character(len=*),           intent(in)  :: shortname
    character(len=*), optional, intent(out) :: longname
    character(len=*), optional, intent(out) :: stdname
    character(len=*), optional, intent(out) :: units

    ! local variables
    integer                    :: i,n
    character(len=CSS)         :: llongname, lstdname, lunits, lshortname  ! local copies
    character(len=*),parameter :: unknown = 'unknown'
    logical                    :: found
    character(len=*),parameter :: subname = '(shr_nuopc_fldList_GetMetadata) '
    ! ----------------------------------------------

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

  subroutine shr_nuopc_fldList_AddFld(flds, stdname, shortname)

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

    ! local variables
    integer :: n,oldsize,id
    logical :: found
    type(shr_nuopc_fldList_entry_type), pointer :: newflds(:)
    character(len=*), parameter :: subname='(shr_nuopc_fldList_AddFld)'
    ! ----------------------------------------------

    if (associated(flds)) then
       oldsize = size(flds)
       found = .false.
       do n= 1,oldsize
          if (trim(stdname) == trim(flds(n)%stdname)) then
             found = .true.
             exit
          end if
       end do
    else
       oldsize = 0
       found = .false.
    end if
    id = oldsize + 1

    ! create new entry if fldname is not in original list

    if (.not. found) then

       ! 1) allocate newfld to be size (one element larger than input flds)
       allocate(newflds(id))

       ! 2) copy flds into first N-1 elements of newflds
       do n = 1,oldsize
          newflds(n)%stdname            = flds(n)%stdname
          newflds(n)%shortname          = flds(n)%shortname
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
    end if

  end subroutine shr_nuopc_fldList_AddFld

  !================================================================================

  subroutine shr_nuopc_fldList_AddMrg(flds, fldname, &
       mrg_from1, mrg_fld1, mrg_type1, mrg_fracname1, &
       mrg_from2, mrg_fld2, mrg_type2, mrg_fracname2, &
       mrg_from3, mrg_fld3, mrg_type3, mrg_fracname3, &
       mrg_from4, mrg_fld4, mrg_type4, mrg_fracname4)

    ! ----------------------------------------------
    ! Determine mrg entry or entries in flds aray
    ! ----------------------------------------------

    type(shr_nuopc_fldList_entry_type) , pointer                :: flds(:)
    character(len=*)                   , intent(in)             :: fldname
    integer                            , intent(in)  , optional :: mrg_from1
    character(len=*)                   , intent(in)  , optional :: mrg_fld1
    character(len=*)                   , intent(in)  , optional :: mrg_type1
    character(len=*)                   , intent(in)  , optional :: mrg_fracname1
    integer                            , intent(in)  , optional :: mrg_from2
    character(len=*)                   , intent(in)  , optional :: mrg_fld2
    character(len=*)                   , intent(in)  , optional :: mrg_type2
    character(len=*)                   , intent(in)  , optional :: mrg_fracname2
    integer                            , intent(in)  , optional :: mrg_from3
    character(len=*)                   , intent(in)  , optional :: mrg_fld3
    character(len=*)                   , intent(in)  , optional :: mrg_type3
    character(len=*)                   , intent(in)  , optional :: mrg_fracname3
    integer                            , intent(in)  , optional :: mrg_from4
    character(len=*)                   , intent(in)  , optional :: mrg_fld4
    character(len=*)                   , intent(in)  , optional :: mrg_type4
    character(len=*)                   , intent(in)  , optional :: mrg_fracname4

    ! local variables
    integer :: n, id
    character(len=*), parameter :: subname='(shr_nuopc_fldList_MrgFld)'
    ! ----------------------------------------------

    id = 0
    do n= 1,size(flds)
       if (trim(fldname) == trim(flds(n)%stdname)) then
          id = n
          exit
       end if
    end do
    if (id == 0) then
       do n = 1,size(flds)
          write(6,*) trim(subname)//' input flds entry is ',trim(flds(n)%stdname)
       end do
       call shr_sys_abort(subname // 'ERROR: fldname '// trim(fldname) // ' not found in input flds')
    end if

    if (present(mrg_from1) .and. present(mrg_fld1) .and. present(mrg_type1)) then
       n = mrg_from1
       flds(id)%merge_fields(n) = mrg_fld1
       flds(id)%merge_types(n) = mrg_type1
       if (present(mrg_fracname1)) then
          flds(id)%merge_fracnames(n) = mrg_fracname1
       end if
    end if
    if (present(mrg_from2) .and. present(mrg_fld2) .and. present(mrg_type2)) then
       n = mrg_from2
       flds(id)%merge_fields(n) = mrg_fld2
       flds(id)%merge_types(n) = mrg_type2
       if (present(mrg_fracname2)) then
          flds(id)%merge_fracnames(n) = mrg_fracname2
       end if
    end if
    if (present(mrg_from3) .and. present(mrg_fld3) .and. present(mrg_type3)) then
       n = mrg_from3
       flds(id)%merge_fields(n) = mrg_fld3
       flds(id)%merge_types(n) = mrg_type3
       if (present(mrg_fracname3)) then
          flds(id)%merge_fracnames(n) = mrg_fracname3
       end if
    end if
    if (present(mrg_from4) .and. present(mrg_fld4) .and. present(mrg_type4)) then
       n = mrg_from4
       flds(id)%merge_fields(n) = mrg_fld4
       flds(id)%merge_types(n) = mrg_type4
       if (present(mrg_fracname4)) then
          flds(id)%merge_fracnames(n) = mrg_fracname4
       end if
    end if

  end subroutine shr_nuopc_fldList_AddMrg

  !================================================================================

  subroutine shr_nuopc_fldList_AddMap(flds, fldname, destcomp, maptype, mapnorm, mapfile)

    ! intput/output variables
    type(shr_nuopc_fldList_entry_type) , intent(inout) :: flds(:)
    character(len=*)                   , intent(in)    :: fldname
    integer                            , intent(in)    :: destcomp
    integer                            , intent(in)    :: maptype
    character(len=*)                   , intent(in)    :: mapnorm
    character(len=*)                   , intent(in)    :: mapfile

    ! local variables
    integer :: id, n
    character(len=*),parameter  :: subname='(shr_nuopc_fldList_AddMap)'
    ! ----------------------------------------------

    id = 0
    do n = 1,size(flds)
       if (trim(fldname) == trim(flds(n)%stdname)) then
          id = n
          exit
       end if
    end do
    if (id == 0) then
       do n = 1,size(flds)
          write(6,*) trim(subname)//' input flds entry is ',trim(flds(n)%stdname)
       end do
       call shr_sys_abort(subname // 'ERROR: fldname '// trim(fldname) // ' not found in input flds')
    end if

    ! Note - default values are already set for the fld entries - so only non-default
    ! values need to be set below
    ! If mapindex is mapfcopy - create a redistribution route handle
    ! If mapfile is idmap - create a redistribution route nhandle
    ! If mapfile is unset then create the mapping route handle at run time

    flds(id)%mapindex(destcomp) = maptype
    flds(id)%mapnorm(destcomp)  = trim(mapnorm)
    flds(id)%mapfile(destcomp)  = trim(mapfile)

    ! overwrite values if appropriate
    if (flds(id)%mapindex(destcomp) == mapfcopy) then
       flds(id)%mapfile(destcomp) = 'unset'
       flds(id)%mapnorm(destcomp) = 'unset'
    else if (trim(flds(id)%mapfile(destcomp)) == 'idmap') then
       flds(id)%mapindex(destcomp) = mapfcopy
       flds(id)%mapnorm(destcomp) = 'unset'
    end if

  end subroutine shr_nuopc_fldList_AddMap

  !================================================================================

  subroutine shr_nuopc_fldList_Realize(state, fldList, flds_scalar_name, flds_scalar_num, &
       grid, mesh, tag, rc)

    use NUOPC             , only : NUOPC_GetStateMemberLists, NUOPC_IsConnected, NUOPC_Realize
    use NUOPC             , only : NUOPC_GetAttribute
    use ESMF              , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF              , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Grid, ESMF_Mesh
    use ESMF              , only : ESMF_StateGet, ESMF_LogFoundError
    use ESMF              , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_FAILURE, ESMF_LOGERR_PASSTHRU
    use ESMF              , only : ESMF_LOGMSG_INFO, ESMF_StateRemove, ESMF_SUCCESS
    use med_constants_mod , only : dbug_flag=>med_constants_dbug_flag

    ! input/output variables
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
           ungriddedUBound=(/flds_scalar_num/), &
           gridToFieldMap=(/2/), &
           rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine shr_nuopc_fldList_Realize

  !================================================================================

  subroutine shr_nuopc_fldList_GetFldInfo_general(fldList, fldindex, stdname, shortname)
    ! ----------------------------------------------
    ! Get field info
    ! ----------------------------------------------
    type(shr_nuopc_fldList_type) , intent(in)  :: fldList
    integer                      , intent(in)  :: fldindex
    character(len=*)             , intent(out) :: stdname
    character(len=*)             , intent(out) :: shortname

    ! local variables
    character(len=*), parameter :: subname='(shr_nuopc_fldList_GetFldInfo_general)'
    ! ----------------------------------------------

    stdname   = fldList%flds(fldindex)%stdname
    shortname = fldList%flds(fldindex)%shortname

  end subroutine shr_nuopc_fldList_GetFldInfo_general

  !================================================================================

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

  !================================================================================

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

    ! input/output variables
    type(shr_nuopc_fldList_type), intent(in)  :: fldList
    ! ----------------------------------------------

    if (associated(fldList%flds)) then
       shr_nuopc_fldList_GetNumFlds = size(fldList%flds)
    else
       shr_nuopc_fldList_GetNumFlds = 0
    end if

  end function shr_nuopc_fldList_GetNumFlds

  !================================================================================

  subroutine shr_nuopc_fldList_GetFldNames(flds, fldnames, rc)

    use ESMF, only : ESMF_LOGMSG_INFO, ESMF_FAILURE, ESMF_SUCCESS, ESMF_LogWrite

    ! input/output variables
    type(shr_nuopc_fldList_entry_type) , pointer     :: flds(:)
    character(len=*)                   , pointer     :: fldnames(:)
    integer, optional                  , intent(out) :: rc

    !local variables
    integer :: n
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (associated(flds) .and. associated(fldnames)) then
       do n = 1,size(flds)
          fldnames(n) = trim(flds(n)%shortname)
       end do
    else
       call ESMF_LogWrite("shr_nuopc_fldList_GetFldNames: ERROR either flds or fldnames have not been allocate ", &
            ESMF_LOGMSG_INFO, rc=rc)
       rc=ESMF_FAILURE
       return
    end if

  end subroutine shr_nuopc_fldList_GetFldNames

  !================================================================================

  subroutine shr_nuopc_fldList_Document_Mapping(logunit, med_coupling_active)

    ! input/output variables
    integer, intent(in)  :: logunit
    logical, intent(in)  :: med_coupling_active(:,:) 

    ! local variables
    integer           :: nsrc,ndst,nf,nm,n
    integer           :: mapindex
    character(len=CS) :: mapnorm
    character(len=CL) :: mapfile
    character(len=CS) :: fldname
    character(len=CS) :: stdname
    character(len=CX) :: merge_fields
    character(len=CX) :: merge_field
    character(len=CS) :: merge_type
    character(len=CS) :: merge_fracname
    character(len=CS) :: string
    character(len=CL) :: mrgstr
    character(len=CL) :: cvalue
    logical           :: init_mrgstr
    character(len=*),parameter :: subname = '(shr_nuopc_fldList_Document_Mapping)'
    !-----------------------------------------------------------

    !---------------------------------------
    ! Document mapping (also add albedo and aoflux) 
    !---------------------------------------

    ! Loop over src components
    do nsrc = 1,ncomps
       ! Loop over all possible destination components for each src component
       do ndst = 1,ncomps
          if (nsrc /= ndst .and. med_coupling_active(nsrc,ndst)) then
             ! Write all the mappings for fields from the src to the destination component
             write(logunit,*)' '
             do n = 1,size(fldListFr(nsrc)%flds)
                mapindex = fldListFr(nsrc)%flds(n)%mapindex(ndst)
                if ( mapindex /= mapunset) then
                   fldname  = trim(fldListFr(nsrc)%flds(n)%stdname)
                   mapnorm  = trim(fldListFr(nsrc)%flds(n)%mapnorm(ndst))
                   mapfile  = trim(fldListFr(nsrc)%flds(n)%mapfile(ndst))

                   if (trim(mapnorm) == 'unset') then
                      cvalue = ' mapping '//trim(compname(nsrc))//'->'//trim(compname(ndst)) //' '//trim(fldname) // &
                           ' via '// trim(mapnames(mapindex))
                   else
                      cvalue = ' mapping '//trim(compname(nsrc))//'->'//trim(compname(ndst)) //' '//trim(fldname) // &
                           ' via '// trim(mapnames(mapindex)) // ' with '// trim(mapnorm) // ' normalization'
                   end if
                   write(logunit,100) trim(cvalue)
                   if (trim(mapfile) /= 'unset' .and. trim(mapfile) /= 'idmap') then
                      cvalue = ' and the mapping file '// trim(mapfile)
                      write(logunit,101) trim(cvalue)
                   end if
                end if
             end do

          end if
       end do
    end do

    ! ocn-> atm mappings for atm/ocn fluxes computed in mediator on the ocn grid
    nsrc = compocn
    ndst = compatm
    if (med_coupling_active(nsrc,ndst)) then
       do n = 1,size(fldListMed_aoflux%flds)
          mapindex = fldlistMed_aoflux%flds(n)%mapindex(ndst)
          if ( mapindex /= mapunset) then
             fldname  = trim(fldlistMed_aoflux%flds(n)%stdname)
             mapnorm  = trim(fldlistMed_aoflux%flds(n)%mapnorm(ndst))
             mapfile  = trim(fldlistMed_aoflux%flds(n)%mapfile(ndst))

             if (trim(mapnorm) == 'unset') then
                cvalue = ' mapping '//trim(compname(nsrc))//'->'//trim(compname(ndst)) //' '//trim(fldname) // &
                     ' via '// trim(mapnames(mapindex))
             else
                cvalue = ' mapping '//trim(compname(nsrc))//'->'//trim(compname(ndst)) //' '//trim(fldname) // &
                     ' via '// trim(mapnames(mapindex)) // ' with '// trim(mapnorm) // ' normalization'
             end if
             write(logunit,100) trim(cvalue)
             if (trim(mapfile) /= 'unset' .and. trim(mapfile) /= 'idmap') then
                cvalue = ' and the mapping file '// trim(mapfile)
                write(logunit,101) trim(cvalue)
             end if
          end if
       end do
    end if

100 format(a)
101 format(3x,a)

  end subroutine shr_nuopc_fldList_Document_Mapping

  !================================================================================

  subroutine shr_nuopc_fldList_Document_Merging(logunit, med_coupling_active)

    !---------------------------------------
    ! Document merging to target destination fields 
    !---------------------------------------

    ! input/output variables
    integer, intent(in)  :: logunit
    logical, intent(in)  :: med_coupling_active(:,:) 

    ! local variables
    integer           :: nsrc,ndst,nf,n
    character(len=CS) :: dst_comp
    character(len=CS) :: dst_field
    character(len=CS) :: src_comp
    character(len=CS) :: src_field
    character(len=CS) :: merge_type
    character(len=CS) :: merge_field
    character(len=CS) :: merge_frac
    character(len=CS) :: prefix
    character(len=CS) :: string
    character(len=CL) :: mrgstr
    logical           :: init_mrgstr
    character(len=*),parameter :: subname = '(shr_nuopc_fldList_Document_Mapping)'
    !-----------------------------------------------------------

    write(logunit,*)

    ! Loop over destination components
    do ndst = 1,ncomps
       dst_comp = trim(compname(ndst))
       prefix = '(merge_to_'//trim(dst_comp)//')'

       ! Loop over all flds in the destination component and determine merging data
       do nf = 1,size(fldListTo(ndst)%flds)
          dst_field = fldListTo(ndst)%flds(nf)%stdname

          ! Loop over all possible source components for destination component field
          mrgstr = ' '
          do nsrc = 1,ncomps

             if (nsrc /= ndst .and. med_coupling_active(nsrc,ndst)) then
                src_comp    = compname(nsrc)
                merge_field = fldListTo(ndst)%flds(nf)%merge_fields(nsrc)
                merge_type  = fldListTo(ndst)%flds(nf)%merge_types(nsrc)
                merge_frac  = fldListTo(ndst)%flds(nf)%merge_fracnames(nsrc)

                if (merge_type == 'merge' .or. merge_type == 'sum_with_weights') then
                   string = trim(merge_frac)//'*'//trim(merge_field)//'('//trim(src_comp)//')'                
                   if (mrgstr == ' ') then
                      mrgstr = trim(prefix)//": "// trim(dst_field) //'('//trim(dst_comp)//')'//' = '//trim(string)
                   else
                      mrgstr = trim(mrgstr) //' + '//trim(string)
                   end if
                else if (merge_type == 'sum') then
                   string = trim(merge_field)//'('//trim(src_comp)//')'
                   if (mrgstr == ' ') then
                      mrgstr = trim(prefix)//": "//trim(dst_field) //'('//trim(dst_comp)//')'//' = '//trim(string)
                   else
                      mrgstr = trim(mrgstr) //' + '//trim(string)
                   end if
                else
                   if (merge_type == 'copy') then
                      mrgstr = trim(prefix)//": " // trim(dst_field) //'('//trim(dst_comp)//')'//' = '// &
                           trim(merge_field)//'('//trim(src_comp)//')'
                   else if (merge_type == 'copy_with_weights') then
                      mrgstr = trim(prefix)//": "// trim(dst_field) //'('//trim(dst_comp)//')'//' = '// &
                           trim(merge_frac)//'*'//trim(merge_field)//'('//trim(src_comp)//')'
                   end if
                end if
             end if

          end do ! end loop over nsrc
          if (mrgstr /= ' ') then
             write(logunit,'(a)') trim(mrgstr)
          end if
       end do ! end loop over nf
       !write(logunit,*)' '
    end do  ! end loop over ndst

  end subroutine shr_nuopc_fldList_Document_Merging

end module esmflds
