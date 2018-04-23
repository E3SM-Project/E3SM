module med_merge_new_mod

  !-----------------------------------------------------------------------------
  ! Performs merges from source field bundles to destination field bundle
  !-----------------------------------------------------------------------------

  use ESMF
  use shr_kind_mod  
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_type
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_GetNumFlds
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_GetFldInfo
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetNameN
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use med_constants_mod

  implicit none 
  private

  integer                       :: dbrc
  integer           , parameter :: dbug_flag   = med_constants_dbug_flag
  real(ESMF_KIND_R8), parameter :: spval_init  = med_constants_spval_init
  real(ESMF_KIND_R8), parameter :: spval       = med_constants_spval
  real(ESMF_KIND_R8), parameter :: czero       = med_constants_czero
  character(*),parameter :: u_FILE_u = &
       __FILE__

  public  :: med_merge_auto
  private :: med_merge

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_merge_auto(FBOut, FBfrac, FBImp, fldListTo, rc)
    type(ESMF_FieldBundle)       , intent(inout) :: FBOut     ! Merged output field bundle
    type(ESMF_FieldBundle)       , intent(in)    :: FBfrac    ! Fraction data for FBOut
    type(ESMF_FieldBundle)       , intent(in)    :: FBImp(:)  ! Array of field bundles each mapping to the FBOut mesh
    type(shr_nuopc_fldList_type) , intent(in)    :: fldListTo ! Information for merging
    integer                      , intent(out)   :: rc

    ! local variables
    integer                :: cnt
    integer                :: n,nf,compsrc
    character(SHR_KIND_CX) :: fldname, stdname
    character(SHR_KIND_CX) :: merge_field
    character(SHR_KIND_CS) :: merge_type
    character(SHR_KIND_CS) :: merge_fracname
    logical,save           :: first_call = .true.
    character(len=*),parameter  :: subname='(med_merge_auto)'
    !---------------------------------------

    ! if (dbug_flag > 5) then
    !    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    ! endif
    rc = ESMF_SUCCESS

    call shr_nuopc_methods_FB_reset(FBOut, value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Want to loop over all of the fields in FBout here - and find the corresponding index in fldListTo(complnd)
    ! for that field name - then call the corresponding merge routine below appropriately

    call ESMF_FieldBundleGet(FBOut, fieldCount=cnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Loop over all fields in field bundle FBOut
    do n = 1,cnt

       ! Get the nth field name in FBexp
       call shr_nuopc_methods_FB_getNameN(FBOut, n, fldname, rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Loop over the field in fldListTo
       do nf = 1,shr_nuopc_fldList_GetNumFlds(fldListTo)

          ! Determine if if there is a match of the fldList field name with the FBOut field name
          call shr_nuopc_fldList_GetFldInfo(fldListTo, nf, stdname)

          if (trim(stdname) == trim(fldname)) then

             ! Loop over all possible source components in the merging arrays returned from the above call
             ! If the merge field name from the source components is not set, then simply go to the next component
             do compsrc = 1,size(FBImp)

                ! First check if the target field bundle to merge exists
                if (ESMF_FieldBundleIsCreated(FBImp(compsrc), rc=rc)) then

                   ! Field bundle exists and field was found -  Now determine the merge information for that field
                   call shr_nuopc_fldList_GetFldInfo(fldListTo, nf, compsrc, merge_field, merge_type, merge_fracname)
                   
                   if (merge_type /= 'unset' .and. merge_field /= 'unset') then
                      if (merge_type == 'copy') then
                         write(6,"(a,i4,a,a,a,a)")'DEBUG: compsrc= ',compsrc, ' merge_type= ',trim(merge_type),' merge_field= ',trim(merge_field) 
                         call med_merge('copy', &
                              FBOut, fldname, &
                              FB=FBImp(compsrc), FBfld=merge_field, &
                              document=first_call, string=subname, rc=rc)
                         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                         
                      else if (merge_type == 'copy_with_weights') then
                         call med_merge('copy', &
                              FBOut, fldname, &
                              FB=FBImp(compsrc), FBFld=merge_field, &
                              FBw=FBfrac, fldw=trim(merge_fracname), &
                              document=first_call, string=subname, rc=rc)
                         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                         
                      else if (merge_type == 'merge') then
                         call med_merge('merge', &
                              FBOut, fldname, &
                              FB=FBImp(compsrc), FBFld=merge_field, &
                              FBw=FBfrac, fldw=trim(merge_fracname), &
                              document=first_call, string=subname, rc=rc)
                         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                      end if
                   end if
                end if
             end do
          end if
       end do
    end do

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    ! if (dbug_flag > 5) then
    !    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    ! endif

  end subroutine med_merge_auto

  !-----------------------------------------------------------------------------

  subroutine med_merge(merge_type, FBout, FBoutfld, FB, FBfld, FBw, fldw, document, string, rc)

    character(len=*)      ,intent(in)          :: merge_type
    type(ESMF_FieldBundle),intent(inout)       :: FBout
    character(len=*)      ,intent(in)          :: FBoutfld
    type(ESMF_FieldBundle),intent(in)          :: FB
    character(len=*)      ,intent(in)          :: FBfld  
    type(ESMF_FieldBundle),intent(in),optional :: FBw
    character(len=*)      ,intent(in),optional :: fldw
    logical               ,intent(in),optional :: document
    character(len=*)      ,intent(in),optional :: string
    integer               ,intent(out)         :: rc

    ! local variables
    real(ESMF_KIND_R8), pointer :: dp1 (:), dp2(:,:)
    real(ESMF_KIND_R8), pointer :: dpf1(:), dpf2(:,:)
    real(ESMF_KIND_R8), pointer :: dpw1(:), dpw2(:,:)
    integer                     :: nf
    integer                     :: lrank
    logical                     :: ldocument
    character(len=128)          :: lstring
    character(len=*),parameter  :: subname='(med_merge)'
    !---------------------------------------

    ! if (dbug_flag > 5) then
    !    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    ! endif
    rc = ESMF_SUCCESS

    ldocument = .false.
    if (present(document)) ldocument=document

    lstring = ""
    if (present(string)) lstring=string

    !-------------------------
    ! Error checks
    !-------------------------

    if ((present(FBw) .and. .not.present(fldw)) .or. (.not.present(FBw) .and. present(fldw))) then
       call ESMF_LogWrite(trim(subname)//": error FBw and fldw both required", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
       rc = ESMF_FAILURE
       return
    end if

    if (present(fldw)) then
       if (trim(fldw) == 'unset') then
          call ESMF_LogWrite(trim(subname)//": error required merge_fracname is not set", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
       end if
    end if

    !-------------------------
    ! Get appropriate field pointers
    !-------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBout, trim(FBoutfld), fldptr1=dp1, fldptr2=dp2, rank=lrank, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (present(FBw) .and. present(fldw)) then
       if (lrank == 1) then
          call shr_nuopc_methods_FB_GetFldPtr(FBw, trim(fldw), fldptr1=dpw1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (lrank == 2) then
          call shr_nuopc_methods_FB_GetFldPtr(FBw, trim(fldw), fldptr2=dpw2, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    endif

    !-------------------------
    ! Loop over all output fields and do the merge
    !-------------------------

    ! Get field pointer to input field used in the merge
    if (lrank == 1) then
       write(6,*)'DEBUG: fbfld = ',trim(fbfld)
       call shr_nuopc_methods_FB_GetFldPtr(FB, trim(FBfld), fldptr1=dpf1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (lrank == 2) then
       call shr_nuopc_methods_FB_GetFldPtr(FB, trim(FBfld), fldptr2=dpf2, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Do one of two types of merges (copy or merge)
    if (trim(merge_type)  == 'copy') then
       if (present(FBw) .and. present(fldw)) then
          if (document) then
             call ESMF_LogWrite(trim(subname)//":"//trim(lstring)//":"//trim(FBoutfld)//&
                  " ="//trim(FBfld), ESMF_LOGMSG_INFO, rc=dbrc)
          end if
          if (lrank == 1) then
             dp1(:) = dpf1(:)*dpw1(:)
          else
             dp2(:,:) = dpf2(:,:)*dpw2(:,:)
          endif
       else
          if (lrank == 1) then
             dp1(:) = dpf1(:)
          else
             dp2(:,:) = dpf2(:,:)
          endif
       end if
    else if (trim(merge_type)  == 'merge') then
       if (present(FBw) .and. present(fldw)) then
          if (document) call ESMF_LogWrite(trim(subname)//":"//trim(lstring)//":"//trim(FBoutfld)//&
               "+="//trim(FBfld)//"*"//trim(fldw), ESMF_LOGMSG_INFO, rc=dbrc)
          if (lrank == 1) then
             dp1(:) = dp1(:) + dpf1(:)*dpw1(:)
          else
             dp2(:,:) = dp2(:,:) + dpf2(:,:)*dpw2(:,:)
          endif
       else
          if (document) call ESMF_LogWrite(trim(subname)//":"//trim(lstring)//":"//trim(FBoutfld)//&
               "+="//trim(FBfld), ESMF_LOGMSG_INFO, rc=dbrc)
          if (lrank == 1) then
             dp1(:) = dp1(:) + dpf1(:)
          else
             dp2(:,:) = dp2(:,:) + dpf2(:,:)
          endif
       endif
    else
       call ESMF_LogWrite(trim(subname)//": merge type "//trim(merge_type)//" not supported", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
       rc = ESMF_FAILURE
       return
    end if

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    ! if (dbug_flag > 5) then
    !    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    ! endif

  end subroutine med_merge

end module med_merge_new_mod
