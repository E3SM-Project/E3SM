module med_merge_mod

  !-----------------------------------------------------------------------------
  ! Performs merges from source field bundles to destination field bundle
  !-----------------------------------------------------------------------------

  use med_constants_mod     , only : R8
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
  use med_constants_mod     , only : spval_init => med_constants_spval_init
  use med_constants_mod     , only : spval => med_constants_spval
  use med_constants_mod     , only : czero => med_constants_czero
  use shr_nuopc_methods_mod , only : ChkErr => shr_nuopc_methods_ChkErr

  implicit none
  private

  public  :: med_merge_auto
  public  :: med_merge_field

  interface med_merge_field ; module procedure &
       med_merge_field_1D, &
       med_merge_field_2D
  end interface

  private :: med_merge_auto_field

  character(*),parameter :: u_FILE_u = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_merge_auto(compout_name, FBOut, FBfrac, FBImp, fldListTo, FBMed1, FBMed2, &
       document, string, mastertask, rc)

    use ESMF                  , only : ESMF_FieldBundle
    use ESMF                  , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogWrite, ESMF_LogMsg_Info
    use ESMF                  , only : ESMF_LogSetError, ESMF_RC_OBJ_NOT_CREATED
    use med_constants_mod     , only : CL, CX, CS
    use shr_string_mod        , only : shr_string_listGetNum
    use shr_string_mod        , only : shr_string_listGetName
    use esmFlds               , only : compmed, compname
    use esmFlds               , only : shr_nuopc_fldList_type
    use esmFlds               , only : shr_nuopc_fldList_GetNumFlds
    use esmFlds               , only : shr_nuopc_fldList_GetFldInfo
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FldChk
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetNameN
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use med_internalstate_mod , only : logunit
    use perf_mod              , only : t_startf, t_stopf

    ! ----------------------------------------------
    ! Auto merge based on fldListTo info
    ! ----------------------------------------------

    ! input/output variables
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
    integer       :: cnt
    integer       :: n,nf,nm,compsrc
    character(CX) :: fldname, stdname
    character(CX) :: merge_fields
    character(CX) :: merge_field
    character(CS) :: merge_type
    character(CS) :: merge_fracname
    character(CL) :: mrgstr   ! temporary string
    logical       :: init_mrgstr
    integer       :: dbrc
    character(len=*),parameter :: subname=' (module_med_merge_mod: med_merge_auto)'
    !---------------------------------------
    call t_startf('MED:'//subname)

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    rc = ESMF_SUCCESS

    call shr_nuopc_methods_FB_reset(FBOut, value=czero, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Want to loop over all of the fields in FBout here - and find the corresponding index in fldListTo(compxxx)
    ! for that field name - then call the corresponding merge routine below appropriately

    call ESMF_FieldBundleGet(FBOut, fieldCount=cnt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (document) then
       if (mastertask) write(logunit,*) ' '
    end if

    ! Loop over all fields in field bundle FBOut
    do n = 1,cnt

       ! Get the nth field name in FBexp
       call shr_nuopc_methods_FB_getNameN(FBOut, n, fldname, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

                ! If merge_field is a colon delimited string then cycle through every field - otherwise by default nm
                ! will only equal 1
                do nm = 1,shr_string_listGetNum(merge_fields)

                   call shr_string_listGetName(merge_fields, nm, merge_field)

                   if (merge_type /= 'unset' .and. merge_field /= 'unset') then

                      write(6,*)'DEBUG: nf, merge_type, merge_fields=', nf, trim(merge_type), trim(merge_fields)

                      ! Document merging if appropriate
                      if (document) then
                         if (merge_type == 'merge' .or. merge_type == 'sum_with_weights') then
                            if (init_mrgstr) then
                               mrgstr = trim(string)//": "// trim(fldname) //'('//trim(compout_name)//')'//' = ' &
                                    // trim(merge_fracname)//'*'//trim(merge_field)//'('//trim(compname(compsrc))//')'
                               init_mrgstr = .false.
                            else
                               mrgstr = trim(mrgstr) //' + ' &
                                    // trim(merge_fracname)//'*'//trim(merge_field)//'('//trim(compname(compsrc))//')'
                            end if
                         else if (merge_type == 'sum') then
                            if (init_mrgstr) then
                               mrgstr = trim(string)//": "// trim(fldname) //'('//trim(compout_name)//')'//' = ' &
                                    //trim(merge_field) //'('//trim(compname(compsrc))//')'
                            else
                               mrgstr = trim(mrgstr) //' + '//trim(merge_field)//'('//trim(compname(compsrc))//')'
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
                            if (.not. ESMF_FieldBundleIsCreated(FBMed1)) then
                               call ESMF_LogSetError(ESMF_RC_OBJ_NOT_CREATED,  &
                                    msg="Field bundle FBMed1 not created.", &
                                    line=__LINE__, file=u_FILE_u, rcToReturn=rc)
                               return
                            endif
                            if (.not. ESMF_FieldBundleIsCreated(FBMed2)) then
                               call ESMF_LogSetError(ESMF_RC_OBJ_NOT_CREATED,  &
                                    msg="Field bundle FBMed2 not created.", &
                                    line=__LINE__, file=u_FILE_u, rcToReturn=rc)
                               return
                            endif
                            if (shr_nuopc_methods_FB_FldChk(FBMed1, trim(merge_field), rc=rc)) then
                               call med_merge_auto_field(trim(merge_type), &
                                    FBOut, fldname, FB=FBMed1, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                               if (ChkErr(rc,__LINE__,u_FILE_u)) return

                            else if (shr_nuopc_methods_FB_FldChk(FBMed2, trim(merge_field), rc=rc)) then
                               call med_merge_auto_field(trim(merge_type), &
                                    FBOut, fldname, FB=FBMed2, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                               if (ChkErr(rc,__LINE__,u_FILE_u)) return

                            else
                               call ESMF_LogWrite(trim(subname)//": ERROR merge_field = "//trim(merge_field)//" not found", &
                                    ESMF_LOGMSG_INFO, rc=rc)
                               rc = ESMF_FAILURE
                               if (ChkErr(rc,__LINE__,u_FILE_u)) return
                            end if

                         elseif (present(FBMed1)) then
                            if (.not. ESMF_FieldBundleIsCreated(FBMed1)) then
                               call ESMF_LogSetError(ESMF_RC_OBJ_NOT_CREATED,  &
                                    msg="Field bundle FBMed1 not created.", &
                                    line=__LINE__, file=u_FILE_u, rcToReturn=rc)
                               return
                            endif
                            if (shr_nuopc_methods_FB_FldChk(FBMed1, trim(merge_field), rc=rc)) then
                               call med_merge_auto_field(trim(merge_type), &
                                    FBOut, fldname, FB=FBMed1, FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                               if (ChkErr(rc,__LINE__,u_FILE_u)) return

                            else
                               call ESMF_LogWrite(trim(subname)//": ERROR merge_field = "//trim(merge_field)//"not found", &
                                    ESMF_LOGMSG_INFO, rc=rc)
                               rc = ESMF_FAILURE
                               if (ChkErr(rc,__LINE__,u_FILE_u)) return
                            end if
                         end if

                      else if (ESMF_FieldBundleIsCreated(FBImp(compsrc), rc=rc)) then

                         if (shr_nuopc_methods_FB_FldChk(FBImp(compsrc), trim(merge_field), rc=rc)) then
                            call med_merge_auto_field(trim(merge_type), &
                                 FBOut, fldname, FB=FBImp(compsrc), FBFld=merge_field, FBw=FBfrac, fldw=trim(merge_fracname), rc=rc)
                            if (ChkErr(rc,__LINE__,u_FILE_u)) return
                         end if
                      end if ! end of single merge

                   end if ! end of check of merge_type and merge_field not unset
                end do ! end of nmerges loop
             end do  ! end of compsrc loop
          end if ! end of check if stdname and fldname are the same
          if (document) then
             if (mastertask) write(logunit,'(a)')trim(mrgstr)
          end if
       end do ! end of loop over fldsListTo
    end do ! end of loop over fields in FBOut

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    call t_stopf('MED:'//subname)

  end subroutine med_merge_auto

  !-----------------------------------------------------------------------------

  subroutine med_merge_auto_field(merge_type, FBout, FBoutfld, FB, FBfld, FBw, fldw, rc)

    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogMsg_Error
    use ESMF                  , only : ESMF_FieldBundle, ESMF_LogWrite, ESMF_LogMsg_Info
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FldChk
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr

    character(len=*)      ,intent(in)    :: merge_type
    type(ESMF_FieldBundle),intent(inout) :: FBout
    character(len=*)      ,intent(in)    :: FBoutfld
    type(ESMF_FieldBundle),intent(in)    :: FB
    character(len=*)      ,intent(in)    :: FBfld
    type(ESMF_FieldBundle),intent(inout) :: FBw
    character(len=*)      ,intent(in)    :: fldw
    integer               ,intent(out)   :: rc

    ! local variables
    real(R8), pointer :: dp1 (:), dp2(:,:)
    real(R8), pointer :: dpf1(:), dpf2(:,:)
    real(R8), pointer :: dpw1(:), dpw2(:,:)
    integer           :: lrank
    integer           :: dbrc
    character(len=*),parameter :: subname=' (med_merge_mod: med_merge)'
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
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (merge_type == 'copy_with_weights' .or. merge_type == 'merge' .or. merge_type == 'sum_with_weights') then
       if (lrank == 1) then
          call shr_nuopc_methods_FB_GetFldPtr(FBw, trim(fldw), fldptr1=dpw1, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (lrank == 2) then
          call shr_nuopc_methods_FB_GetFldPtr(FBw, trim(fldw), fldptr2=dpw2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    endif

    !-------------------------
    ! Loop over all output fields and do the merge
    !-------------------------

    ! Get field pointer to input field used in the merge
    if (lrank == 1) then
       call shr_nuopc_methods_FB_GetFldPtr(FB, trim(FBfld), fldptr1=dpf1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (lrank == 2) then
       call shr_nuopc_methods_FB_GetFldPtr(FB, trim(FBfld), fldptr2=dpf2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
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
    else if (trim(merge_type) == 'sum') then
       if (lrank == 1) then
          dp1(:) = dp1(:) + dpf1(:)
       else
          dp2(:,:) = dp2(:,:) + dpf2(:,:)
       endif
    else if (trim(merge_type) == 'sum_with_weights') then
       if (lrank == 1) then
          dp1(:) = dp1(:) + dpf1(:)*dpw1(:)
       else
          dp2(:,:) = dp2(:,:) + dpf2(:,:)*dpw2(:,:)
       endif
    else
       call ESMF_LogWrite(trim(subname)//": merge type "//trim(merge_type)//" not supported", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
       rc = ESMF_FAILURE
       return
    end if

  end subroutine med_merge_auto_field

  !-----------------------------------------------------------------------------

  subroutine med_merge_field_1D(FBout, fnameout, &
                                FBinA, fnameA, wgtA, &
                                FBinB, fnameB, wgtB, &
                                FBinC, fnameC, wgtC, &
                                FBinD, fnameD, wgtD, &
                                FBinE, fnameE, wgtE, rc)

    use ESMF                  , only : ESMF_FieldBundle, ESMF_LogWrite
    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGMSG_ERROR
    use ESMF                  , only : ESMF_LOGMSG_WARNING, ESMF_LOGMSG_INFO
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FieldPtr_Compare
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FldChk

    ! ----------------------------------------------
    ! Supports up to a five way merge
    ! ----------------------------------------------

    ! input/output variabes
    type(ESMF_FieldBundle) , intent(inout)                 :: FBout
    character(len=*)       , intent(in)                    :: fnameout
    type(ESMF_FieldBundle) , intent(in)                    :: FBinA
    character(len=*)       , intent(in)                    :: fnameA
    real(R8)               , intent(in), pointer           :: wgtA(:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinB
    character(len=*)       , intent(in), optional          :: fnameB
    real(R8)               , intent(in), optional, pointer :: wgtB(:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinC
    character(len=*)       , intent(in), optional          :: fnameC
    real(R8)               , intent(in), optional, pointer :: wgtC(:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinD
    character(len=*)       , intent(in), optional          :: fnameD
    real(R8)               , intent(in), optional, pointer :: wgtD(:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinE
    character(len=*)       , intent(in), optional          :: fnameE
    real(R8)               , intent(in), optional, pointer :: wgtE(:)
    integer                , intent(out)                   :: rc

    ! local variables
    real(R8), pointer          :: dataOut(:)
    real(R8), pointer          :: dataPtr(:)
    real(R8), pointer          :: wgt(:)
    integer                    :: lb1,ub1,i,j,n
    logical                    :: wgtfound, FBinfound
    integer                    :: dbrc
    character(len=*),parameter :: subname='(med_merge_fieldo_1d)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc=ESMF_SUCCESS

    ! check each field has a fieldname passed in
    if ((present(FBinB) .and. .not.present(fnameB)) .or. &
        (present(FBinC) .and. .not.present(fnameC)) .or. &
        (present(FBinD) .and. .not.present(fnameD)) .or. &
        (present(FBinE) .and. .not.present(fnameE))) then

       call ESMF_LogWrite(trim(subname)//": ERROR fname not present with FBin", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
       rc = ESMF_FAILURE
       return
    endif

    if (.not. shr_nuopc_methods_FB_FldChk(FBout, trim(fnameout), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": WARNING field not in FBout, skipping merge "//trim(fnameout), &
            ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
       return
    endif

    call shr_nuopc_methods_FB_GetFldPtr(FBout, trim(fnameout), fldptr1=dataOut, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    lb1 = lbound(dataOut,1)
    ub1 = ubound(dataOut,1)
    allocate(wgt(lb1:ub1))

    dataOut = czero

    ! check that each field passed in actually exists, if not DO NOT do any merge
    FBinfound = .true.
    if (present(FBinB)) then
       if (.not. shr_nuopc_methods_FB_FldChk(FBinB, trim(fnameB), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinC)) then
       if (.not. shr_nuopc_methods_FB_FldChk(FBinC, trim(fnameC), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinD)) then
       if (.not. shr_nuopc_methods_FB_FldChk(FBinD, trim(fnameD), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinE)) then
       if (.not. shr_nuopc_methods_FB_FldChk(FBinE, trim(fnameE), rc=rc)) FBinfound = .false.
    endif
    if (.not. FBinfound) then
       call ESMF_LogWrite(trim(subname)//": WARNING field not found in FBin, skipping merge "//trim(fnameout), &
            ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
       return
    endif

    ! n=1,5 represents adding A to E inputs if they exist
    do n = 1,5
       FBinfound = .false.
       wgtfound = .false.

       if (n == 1) then
          FBinfound = .true.
          call shr_nuopc_methods_FB_GetFldPtr(FBinA, trim(fnameA), fldptr1=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          wgtfound = .true.
          wgt => wgtA

       elseif (n == 2 .and. present(FBinB)) then
          FBinfound = .true.
          call shr_nuopc_methods_FB_GetFldPtr(FBinB, trim(fnameB), fldptr1=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtB)) then
             wgtfound = .true.
             wgt => wgtB
          endif

       elseif (n == 3 .and. present(FBinC)) then
          FBinfound = .true.
          call shr_nuopc_methods_FB_GetFldPtr(FBinC, trim(fnameC), fldptr1=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtC)) then
             wgtfound = .true.
             wgt => wgtC
          endif

       elseif (n == 4 .and. present(FBinD)) then
          FBinfound = .true.
          call shr_nuopc_methods_FB_GetFldPtr(FBinD, trim(fnameD), fldptr1=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtD)) then
             wgtfound = .true.
             wgt => wgtD
          endif

       elseif (n == 5 .and. present(FBinE)) then
          FBinfound = .true.
          call shr_nuopc_methods_FB_GetFldPtr(FBinE, trim(fnameE), fldptr1=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtE)) then
             wgtfound = .true.
             wgt => wgtE
          endif

       endif

       if (FBinfound) then
          if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtr, dataOut, subname, rc)) then
             call ESMF_LogWrite(trim(subname)//": ERROR FBin wrong size", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
             return
          endif

          if (wgtfound) then
             if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtr, wgt, subname, rc)) then
                call ESMF_LogWrite(trim(subname)//": ERROR wgt wrong size", &
                     ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
                rc = ESMF_FAILURE
                return
             endif
             do i = lb1,ub1
                dataOut(i) = dataOut(i) + dataPtr(i) * wgt(i)
             enddo
          else
             do i = lb1,ub1
                dataOut(i) = dataOut(i) + dataPtr(i)
             enddo
          endif  ! wgtfound

       endif  ! FBin found
    enddo  ! n

    if (dbug_flag > 10) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_merge_field_1D

  !-----------------------------------------------------------------------------

  subroutine med_merge_field_2D(FBout, fnameout,     &
                                FBinA, fnameA, wgtA, &
                                FBinB, fnameB, wgtB, &
                                FBinC, fnameC, wgtC, &
                                FBinD, fnameD, wgtD, &
                                FBinE, fnameE, wgtE, rc)

    use ESMF                  , only : ESMF_FieldBundle, ESMF_LogWrite
    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGMSG_ERROR
    use ESMF                  , only : ESMF_LOGMSG_WARNING, ESMF_LOGMSG_INFO
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FieldPtr_Compare
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FldChk

    ! ----------------------------------------------
    ! Supports up to a five way merge
    ! ----------------------------------------------

    ! input/output arguments
    type(ESMF_FieldBundle) , intent(inout)                 :: FBout
    character(len=*)       , intent(in)                    :: fnameout
    type(ESMF_FieldBundle) , intent(in)                    :: FBinA
    character(len=*)       , intent(in)                    :: fnameA
    real(R8)               , intent(in), pointer           :: wgtA(:,:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinB
    character(len=*)       , intent(in), optional          :: fnameB
    real(R8)               , intent(in), optional, pointer :: wgtB(:,:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinC
    character(len=*)       , intent(in), optional          :: fnameC
    real(R8)               , intent(in), optional, pointer :: wgtC(:,:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinD
    character(len=*)       , intent(in), optional          :: fnameD
    real(R8)               , intent(in), optional, pointer :: wgtD(:,:)
    type(ESMF_FieldBundle) , intent(in), optional          :: FBinE
    character(len=*)       , intent(in), optional          :: fnameE
    real(R8)               , intent(in), optional, pointer :: wgtE(:,:)
    integer                , intent(out)                   :: rc

    ! local variables
    real(R8), pointer          :: dataOut(:,:)
    real(R8), pointer          :: dataPtr(:,:)
    real(R8), pointer          :: wgt(:,:)
    integer                    :: lb1,ub1,lb2,ub2,i,j,n
    logical                    :: wgtfound, FBinfound
    integer                    :: dbrc
    character(len=*),parameter :: subname='(med_merge_field_2d)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc=ESMF_SUCCESS

    if (.not. shr_nuopc_methods_FB_FldChk(FBout, trim(fnameout), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": WARNING field not in FBout, skipping merge "//&
            trim(fnameout), ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
       return
    endif

    call shr_nuopc_methods_FB_GetFldPtr(FBout, trim(fnameout), fldptr2=dataOut, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    lb1 = lbound(dataOut,1)
    ub1 = ubound(dataOut,1)
    lb2 = lbound(dataOut,2)
    ub2 = ubound(dataOut,2)
    allocate(wgt(lb1:ub1,lb2:ub2))

    dataOut = czero

    ! check each field has a fieldname passed in
    if ((present(FBinB) .and. .not.present(fnameB)) .or. &
        (present(FBinC) .and. .not.present(fnameC)) .or. &
        (present(FBinD) .and. .not.present(fnameD)) .or. &
        (present(FBinE) .and. .not.present(fnameE))) then
       call ESMF_LogWrite(trim(subname)//": ERROR fname not present with FBin", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
       rc = ESMF_FAILURE
       return
    endif

    ! check that each field passed in actually exists, if not DO NOT do any merge
    FBinfound = .true.
    if (present(FBinB)) then
       if (.not. shr_nuopc_methods_FB_FldChk(FBinB, trim(fnameB), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinC)) then
       if (.not. shr_nuopc_methods_FB_FldChk(FBinC, trim(fnameC), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinD)) then
       if (.not. shr_nuopc_methods_FB_FldChk(FBinD, trim(fnameD), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinE)) then
       if (.not. shr_nuopc_methods_FB_FldChk(FBinE, trim(fnameE), rc=rc)) FBinfound = .false.
    endif
    if (.not. FBinfound) then
       call ESMF_LogWrite(trim(subname)//": WARNING field not found in FBin, skipping merge "//trim(fnameout), &
            ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
       return
    endif

    ! n=1,5 represents adding A to E inputs if they exist
    do n = 1,5
       FBinfound = .false.
       wgtfound = .false.

       if (n == 1) then
          FBinfound = .true.
          call shr_nuopc_methods_FB_GetFldPtr(FBinA, trim(fnameA), fldptr2=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          wgtfound = .true.
          wgt => wgtA

       elseif (n == 2 .and. present(FBinB)) then
          FBinfound = .true.
          call shr_nuopc_methods_FB_GetFldPtr(FBinB, trim(fnameB), fldptr2=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtB)) then
             wgtfound = .true.
             wgt => wgtB
          endif

       elseif (n == 3 .and. present(FBinC)) then
          FBinfound = .true.
          call shr_nuopc_methods_FB_GetFldPtr(FBinC, trim(fnameC), fldptr2=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtC)) then
             wgtfound = .true.
             wgt => wgtC
          endif

       elseif (n == 4 .and. present(FBinD)) then
          FBinfound = .true.
          call shr_nuopc_methods_FB_GetFldPtr(FBinD, trim(fnameD), fldptr2=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtD)) then
             wgtfound = .true.
             wgt => wgtD
          endif

       elseif (n == 5 .and. present(FBinE)) then
          FBinfound = .true.
          call shr_nuopc_methods_FB_GetFldPtr(FBinE, trim(fnameE), fldptr2=dataPtr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(wgtE)) then
             wgtfound = .true.
             wgt => wgtE
          endif

       endif

       if (FBinfound) then
          if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtr, dataOut, subname, rc)) then
             call ESMF_LogWrite(trim(subname)//": ERROR FBin wrong size", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
             return
          endif

          if (wgtfound) then
             if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtr, wgt, subname, rc)) then
                call ESMF_LogWrite(trim(subname)//": ERROR wgt wrong size", &
                     ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
                rc = ESMF_FAILURE
                return
             endif
             do j = lb2,ub2
                do i = lb1,ub1
                   dataOut(i,j) = dataOut(i,j) + dataPtr(i,j) * wgt(i,j)
                enddo
             enddo
          else
             do j = lb2,ub2
                do i = lb1,ub1
                   dataOut(i,j) = dataOut(i,j) + dataPtr(i,j)
                enddo
             enddo
          endif  ! wgtfound

       endif  ! FBin found
    enddo  ! n

    if (dbug_flag > 10) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_merge_field_2D

end module med_merge_mod
