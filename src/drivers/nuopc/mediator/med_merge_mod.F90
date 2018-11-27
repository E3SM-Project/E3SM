module med_merge_mod

  !-----------------------------------------------------------------------------
  ! Performs merges from source field bundles to destination field bundle
  !-----------------------------------------------------------------------------

  use med_constants_mod, only : dbug_flag => med_constants_dbug_flag
  use med_constants_mod, only : spval_init => med_constants_spval_init
  use med_constants_mod, only : spval => med_constants_spval
  use med_constants_mod, only : czero => med_constants_czero

  implicit none
  private

  character(*),parameter :: u_FILE_u = &
       __FILE__

  public  :: med_merge_auto
  private :: med_merge

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_merge_auto(compout_name, FBOut, FBfrac, FBImp, fldListTo, FBMed1, FBMed2, &
       document, string, mastertask, rc)

    use ESMF                  , only : ESMF_FieldBundle
    use ESMF                  , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogWrite, ESMF_LogMsg_Info
    use med_constants_mod     , only : CL, CX, CS
    use shr_string_mod        , only : shr_string_listGetNum
    use shr_string_mod        , only : shr_string_listGetName
    use esmFlds               , only : compmed, compname
    use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_type
    use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_GetNumFlds
    use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_GetFldInfo
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FldChk
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetNameN
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use med_internalstate_mod , only : logunit
    use perf_mod              , only : t_startf, t_stopf
    character(len=*)             , intent(in)            :: compout_name ! component name for FBOut
    type(ESMF_FieldBundle)       , intent(inout)         :: FBOut        ! Merged output field bundle
    type(ESMF_FieldBundle)       , intent(inout)         :: FBfrac       ! Fraction data for FBOut
    type(ESMF_FieldBundle)       , intent(in)            :: FBImp(:)     ! Array of field bundles each mapping to the FBOut mesh
    type(shr_nuopc_fldList_type) , intent(in)            :: fldListTo    ! Information for merging
    type(ESMF_FieldBundle)       , intent(in) , optional :: FBMed1       ! mediator field bundle
    type(ESMF_FieldBundle)       , intent(in) , optional :: FBMed2       ! mediator field bundle
    logical                      , intent(in)            :: document
    character(len=*)             , intent(in)            :: string
    logical                      , intent(in)            :: mastertask
    integer                      , intent(out)           :: rc

    ! local variables
    integer                :: cnt
    integer                :: n,nf,nm,compsrc
    character(CX) :: fldname, stdname
    character(CX) :: merge_fields
    character(CX) :: merge_field
    character(CS) :: merge_type
    character(CS) :: merge_fracname
    character(CL) :: mrgstr   ! temporary string
    logical                :: init_mrgstr
    character(len=*),parameter  :: subname='(med_merge_auto)'
    integer                       :: dbrc
    !---------------------------------------
    call t_startf('MED:'//subname)

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
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

       init_mrgstr = .true.

       ! Loop over the field in fldListTo
       do nf = 1,shr_nuopc_fldList_GetNumFlds(fldListTo)

          ! Determine if if there is a match of the fldList field name with the FBOut field name
          call shr_nuopc_fldList_GetFldInfo(fldListTo, nf, stdname)

          if (trim(stdname) == trim(fldname)) then

             ! Loop over all possible source components in the merging arrays returned from the above call
             ! If the merge field name from the source components is not set, then simply go to the next component
             do compsrc = 1,size(FBImp)

                ! Determine the merge information for the import field
                call shr_nuopc_fldList_GetFldInfo(fldListTo, nf, compsrc, merge_fields, merge_type, merge_fracname)

                do nm = 1,shr_string_listGetNum(merge_fields)
                   call shr_string_listGetName(merge_fields, nm, merge_field)
                   if (merge_type /= 'unset' .and. merge_field /= 'unset') then

                      ! Document merging if appropriate
                      if (document) then
                         if (merge_type == 'merge' .or. merge_type == 'accumulate') then
                            if (init_mrgstr) then
                               mrgstr = trim(string)//": "// trim(fldname) //'('//trim(compout_name)//')'//' = ' &
                                    // trim(merge_fracname)//'*'//trim(merge_field)//'('//trim(compname(compsrc))//')'
                               init_mrgstr = .false.
                            else
                               mrgstr = trim(mrgstr) //' + &
                                    '// trim(merge_fracname)//'*'//trim(merge_field)//'('//trim(compname(compsrc))//')'
                            end if
                         else
                            if (merge_type == 'copy') then
                               mrgstr = trim(string)//": " // trim(fldname) //'('//trim(compout_name)//')'//' = ' &
                                    //trim(merge_field) //'('//trim(compname(compsrc))//')'
                            else if (merge_type == 'copy_with_weights') then
                               mrgstr = trim(string)//": "// trim(fldname) //'('//trim(compout_name)//')'//' = ' &
                                    //trim(merge_fracname)//'*'//trim(merge_field)//'('//trim(compname(compsrc))//')'
                            end if
                         end if
                      end if

                      ! Perform merge
                      if (compsrc == compmed) then

                         if (present(FBMed1) .and. present(FBMed2)) then
                            if (shr_nuopc_methods_FB_FldChk(FBMed1, trim(merge_field), rc=rc)) then
                               call med_merge(trim(merge_type), &
                                    FBOut, fldname, FB=FBMed1, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                            else if (shr_nuopc_methods_FB_FldChk(FBMed2, trim(merge_field), rc=rc)) then
                               call med_merge(trim(merge_type), &
                                    FBOut, fldname, FB=FBMed2, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                            else
                               call ESMF_LogWrite(trim(subname)//": ERROR merge_field = "//trim(merge_field)//"not found", &
                                    ESMF_LOGMSG_INFO, rc=rc)
                               rc = ESMF_FAILURE
                               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                            end if
                         elseif (present(FBMed1)) then
                            if (shr_nuopc_methods_FB_FldChk(FBMed1, trim(merge_field), rc=rc)) then
                               call med_merge(trim(merge_type), &
                                    FBOut, fldname, FB=FBMed1, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                            else
                               call ESMF_LogWrite(trim(subname)//": ERROR merge_field = "//trim(merge_field)//"not found", &
                                    ESMF_LOGMSG_INFO, rc=rc)
                               rc = ESMF_FAILURE
                               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                            end if
                         end if

                      else if (ESMF_FieldBundleIsCreated(FBImp(compsrc), rc=rc)) then

                         if (shr_nuopc_methods_FB_FldChk(FBImp(compsrc), trim(merge_field), rc=rc)) then
                            call med_merge(trim(merge_type), &
                                 FBOut, fldname, FB=FBImp(compsrc), FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                         end if
                      end if ! end of single merge
                   end if ! end of check of merge_type and merge_field not unset
                end do ! end of nmerges loop
             end do  ! end of compsrc loop
             if (document) then
                if (mastertask) write(logunit,'(a)')trim(mrgstr)
             end if
          end if ! end of check if stdname and fldname are the same
       end do ! end of loop over fldsListTo
    end do ! end of loop over fields in FBOut

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    call t_stopf('MED:'//subname)

  end subroutine med_merge_auto

  !-----------------------------------------------------------------------------

  subroutine med_merge(merge_type, FBout, FBoutfld, FB, FBfld, FBw, fldw, rc)
    use ESMF, only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogMsg_Error
    use ESMF, only : ESMF_FieldBundle, ESMF_LogWrite, ESMF_LogMsg_Info
    use med_constants_mod, only : R8
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FldChk
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use med_internalstate_mod , only : logunit

    character(len=*)      ,intent(in)          :: merge_type
    type(ESMF_FieldBundle),intent(inout)       :: FBout
    character(len=*)      ,intent(in)          :: FBoutfld
    type(ESMF_FieldBundle),intent(in)          :: FB
    character(len=*)      ,intent(in)          :: FBfld
    type(ESMF_FieldBundle),intent(inout)          :: FBw !DEBUG - change back to in
    character(len=*)      ,intent(in)          :: fldw
    integer               ,intent(out)         :: rc

    ! local variables
    real(R8), pointer :: dp1 (:), dp2(:,:)
    real(R8), pointer :: dpf1(:), dpf2(:,:)
    real(R8), pointer :: dpw1(:), dpw2(:,:)
    integer                     :: lrank
    character(len=*),parameter  :: subname='(med_merge)'
    integer                       :: dbrc
    !---------------------------------------

    rc = ESMF_SUCCESS

    !-------------------------
    ! Error checks
    !-------------------------

    if (merge_type == 'copy_with_weights' .or. merge_type == 'merge') then
       if (trim(fldw) == 'unset') then
          call ESMF_LogWrite(trim(subname)//": error required merge_fracname is not set", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
       end if
       if (.not. shr_nuopc_methods_FB_FldChk(FBw, trim(fldw), rc=rc)) then
          call ESMF_LogWrite(trim(subname)//": error "//trim(fldw)//"is not in FBw", &
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

    if (merge_type == 'copy_with_weights' .or. merge_type == 'merge') then
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
       call shr_nuopc_methods_FB_GetFldPtr(FB, trim(FBfld), fldptr1=dpf1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (lrank == 2) then
       call shr_nuopc_methods_FB_GetFldPtr(FB, trim(FBfld), fldptr2=dpf2, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Do one of two types of merges (copy or merge)
    if (trim(merge_type)  == 'copy') then
       if (lrank == 1) then
          dp1(:) = dpf1(:)
       else
          dp2(:,:) = dpf2(:,:)
       endif
    else if (trim(merge_type)  == 'copy_with_weights') then
       if (lrank == 1) then
          dp1(:) = dpf1(:)*dpw1(:)
       else
          dp2(:,:) = dpf2(:,:)*dpw2(:,:)
       endif
    else if (trim(merge_type)  == 'merge') then
       if (lrank == 1) then
          dp1(:) = dp1(:) + dpf1(:)*dpw1(:)
       else
          dp2(:,:) = dp2(:,:) + dpf2(:,:)*dpw2(:,:)
       endif
    else if (trim(merge_type) == 'accumulate') then
       if (lrank == 1) then
          dp1(:) = dp1(:) + dpf1(:)
       else
          dp2(:,:) = dp2(:,:) + dpf2(:,:)
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

  end subroutine med_merge

end module med_merge_mod
