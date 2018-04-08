module med_merge_mod

  !-----------------------------------------------------------------------------
  ! Mediator Component.
  ! This mediator operates on two timescales and keeps two internal Clocks to
  ! do so.
  !-----------------------------------------------------------------------------

  use ESMF
  use shr_nuopc_methods_mod
  use med_internalstate_mod
  use med_constants_mod

  implicit none

  private

  integer            :: dbrc
  integer           , parameter :: dbug_flag   = med_constants_dbug_flag
  logical           , parameter :: statewrite_flag = med_constants_statewrite_flag
  real(ESMF_KIND_R8), parameter :: spval_init  = med_constants_spval_init
  real(ESMF_KIND_R8), parameter :: spval       = med_constants_spval
  real(ESMF_KIND_R8), parameter :: czero       = med_constants_czero
  character(len=ESMF_MAXSTR) :: tmpstr
  character(*),parameter :: u_FILE_u = &
    __FILE__

  public med_merge_auto

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

    subroutine med_merge_auto(FBout, &
                              FB1, FB1w, fldw1, &
                              FB2, FB2w, fldw2, &
                              FB3, FB3w, fldw3, &
                              FB4, FB4w, fldw4, &
                              FB5, FB5w, fldw5, &
                              FB6, FB6w, fldw6, &
                              document, string, rc)

    type(ESMF_FieldBundle),intent(inout)       :: FBout
    type(ESMF_FieldBundle),intent(in),optional :: FB1  , FB2  , FB3  , FB4  , FB5  , FB6
    type(ESMF_FieldBundle),intent(in),optional :: FB1w , FB2w , FB3w , FB4w , FB5w , FB6w
    character(len=*)      ,intent(in),optional :: fldw1, fldw2, fldw3, fldw4, fldw5, fldw6
    logical               ,intent(in),optional :: document
    character(len=*)      ,intent(in),optional :: string
    integer               ,intent(out)         :: rc

    ! This subroutine initializes the fractions

    ! local variables
    integer                     :: cnt, n
    logical                     :: ldocument
    character(len=128)          :: lstring
    character(len=*),parameter  :: subname='(med_merge_auto)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    !---------------------------------------

    ldocument = .false.
    if (present(document)) ldocument=document

    lstring = ""
    if (present(string)) lstring=trim(string)

    call ESMF_FieldBundleGet(FBout, fieldCount=cnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,cnt

       if (present(FB1)) then
       if (ESMF_FieldBundleIsCreated(FB1,rc=rc)) then
         if (present(FB1w) .and. present(fldw1)) then
           if (.not.ESMF_FieldBundleIsCreated(FB1w,rc=rc)) then
              call ESMF_LogWrite(trim(subname)//": error FB1w not created", &
                   ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
           else
             call med_merge_fbx(FBout, n, FB1, FB1w, fldw1, document=ldocument, string=trim(lstring), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
           endif
         elseif ((present(FB1w) .and. .not.present(fldw1)) .or. &
                 (.not.present(FB1w) .and. present(fldw1))) then
            call ESMF_LogWrite(trim(subname)//": error FB1w and fldw1 both required", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
           rc = ESMF_FAILURE
           return
         else
           call med_merge_fbx(FBout, n, FB1, document=ldocument, string=trim(lstring), rc=rc)
           if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         endif
       endif
       endif

       if (present(FB2)) then
       if (ESMF_FieldBundleIsCreated(FB2,rc=rc)) then
         if (present(FB2w) .and. present(fldw2)) then
           if (.not.ESMF_FieldBundleIsCreated(FB2w,rc=rc)) then
              call ESMF_LogWrite(trim(subname)//": error FB2w not created", &
                   ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
           else
             call med_merge_fbx(FBout, n, FB2, FB2w, fldw2, document=ldocument, string=trim(lstring), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
           endif
         elseif ((present(FB2w) .and. .not.present(fldw2)) .or. &
                 (.not.present(FB2w) .and. present(fldw2))) then
            call ESMF_LogWrite(trim(subname)//": error FB2w and fldw2 both required", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
           rc = ESMF_FAILURE
           return
         else
           call med_merge_fbx(FBout, n, FB2, document=ldocument, string=trim(lstring), rc=rc)
           if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         endif
       endif
       endif

       if (present(FB3)) then
       if (ESMF_FieldBundleIsCreated(FB3,rc=rc)) then
         if (present(FB3w) .and. present(fldw3)) then
           if (.not.ESMF_FieldBundleIsCreated(FB3w,rc=rc)) then
              call ESMF_LogWrite(trim(subname)//": error FB3w not created", &
                   ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
           else
             call med_merge_fbx(FBout, n, FB3, FB3w, fldw3, document=ldocument, string=trim(lstring), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
           endif
         elseif ((present(FB3w) .and. .not.present(fldw3)) .or. &
                 (.not.present(FB3w) .and. present(fldw3))) then
            call ESMF_LogWrite(trim(subname)//": error FB3w and fldw3 both required", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
           rc = ESMF_FAILURE
           return
         else
           call med_merge_fbx(FBout, n, FB3, document=ldocument, string=trim(lstring), rc=rc)
           if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         endif
       endif
       endif

       if (present(FB4)) then
       if (ESMF_FieldBundleIsCreated(FB4,rc=rc)) then
         if (present(FB4w) .and. present(fldw4)) then
           if (.not.ESMF_FieldBundleIsCreated(FB4w,rc=rc)) then
              call ESMF_LogWrite(trim(subname)//": error FB4w not created", &
                   ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
           else
             call med_merge_fbx(FBout, n, FB4, FB4w, fldw4, document=ldocument, string=trim(lstring), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
           endif
         elseif ((present(FB4w) .and. .not.present(fldw4)) .or. &
                 (.not.present(FB4w) .and. present(fldw4))) then
            call ESMF_LogWrite(trim(subname)//": error FB4w and fldw4 both required", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
           rc = ESMF_FAILURE
           return
         else
           call med_merge_fbx(FBout, n, FB4, document=ldocument, string=trim(lstring), rc=rc)
           if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         endif
       endif
       endif

       if (present(FB5)) then
       if (ESMF_FieldBundleIsCreated(FB5,rc=rc)) then
         if (present(FB5w) .and. present(fldw5)) then
           if (.not.ESMF_FieldBundleIsCreated(FB5w,rc=rc)) then
              call ESMF_LogWrite(trim(subname)//": error FB5w not created", &
                   ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
           else
             call med_merge_fbx(FBout, n, FB5, FB5w, fldw5, document=ldocument, string=trim(lstring), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
           endif
         elseif ((present(FB5w) .and. .not.present(fldw5)) .or. &
                 (.not.present(FB5w) .and. present(fldw5))) then
            call ESMF_LogWrite(trim(subname)//": error FB5w and fldw5 both required", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
           rc = ESMF_FAILURE
           return
         else
           call med_merge_fbx(FBout, n, FB5, document=ldocument, string=trim(lstring), rc=rc)
           if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         endif
       endif
       endif

       if (present(FB6)) then
       if (ESMF_FieldBundleIsCreated(FB6,rc=rc)) then
         if (present(FB6w) .and. present(fldw6)) then
           if (.not.ESMF_FieldBundleIsCreated(FB6w,rc=rc)) then
              call ESMF_LogWrite(trim(subname)//": error FB6w not created", &
                   ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
           else
             call med_merge_fbx(FBout, n, FB6, FB6w, fldw6, document=ldocument, string=trim(lstring), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
           endif
         elseif ((present(FB6w) .and. .not.present(fldw6)) .or. &
                 (.not.present(FB6w) .and. present(fldw6))) then
            call ESMF_LogWrite(trim(subname)//": error FB6w and fldw6 both required", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
           rc = ESMF_FAILURE
           return
         else
           call med_merge_fbx(FBout, n, FB6, document=ldocument, string=trim(lstring), rc=rc)
           if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         endif
       endif
       endif

    enddo

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_merge_auto

  !-----------------------------------------------------------------------------

  subroutine med_merge_fbx(FBout, n, FB, FBw, fldw, document, string, rc)
    type(ESMF_FieldBundle),intent(inout) :: FBout
    integer               ,intent(in)    :: n
    type(ESMF_FieldBundle),intent(in)    :: FB
    type(ESMF_FieldBundle),intent(in),optional :: FBw
    character(len=*)      ,intent(in),optional :: fldw
    logical               ,intent(in),optional :: document
    character(len=*)      ,intent(in),optional :: string
    integer               ,intent(out)         :: rc

    ! This subroutine merges fields from FB into the nth field of FBout.
    ! If FBw and fldw are present, the the FB field is weighted (multiplied) by that field.
    ! document and string are optional arguments related to diagnostics.

    ! A field is expected to have a naming convention like "char1"_"char2"
    ! where the char1 before the underscore is the prefix indicative of State or Flux
    ! and the component associated with the field and char2 is the field name.
    ! Some examples are So_t, Sx_t, Fioi_taux, Faxx_taux.
    ! S is a state, F is a flux.  An x in the prefix indicates it should be merged.
    ! If there is no x in the prefix, then it's a field that is filled via copy.

    ! Here are the merging rules
    ! - The field is copied with no weighting if
    !   - the full fieldname (char1_char2) in FB exactly matches FBout
    ! - The field is merged with weighting (if it exists) if
    !   - there is an "x" in the FBout prefix
    !   - the fieldname (char2) in FB matches the fieldname in FBout
    !   - if the FBout prefix has an F, then either FBpre has an x in it or there is
    !     another character, not x, that matches between the FBout prefix and the FB prefix.

    ! Some examples
    !  FBout = So_t, FB=So_t, copied
    !  FBout = Sx_t, FB=So_t, merged
    !  FBout = So_t, FB=Si_t, skipped
    !  FBout = Fioi_salt, FB=Fioi_salt, copied (identical)
    !  FBout = Foxx_salt, FB=Foxx_salt, copied (identical)
    !  FBout = Foxx_salt, FB=Faxa_salt, merged ("x" in Faxa allows merge)
    !  FBout = Foxx_salt, FB=Fioi_salt, merged ("o" matches in Foxx, Fioi)
    !  FBout = Fioi_salt, FB=Faii_salt, skipped (no "x" in Fioi)
    !  FBout = Foxx_salt, FB=Faii_salt, skipped ("o" in Foxx is not present in Faii)

    ! local variables
    real(ESMF_KIND_R8), pointer :: dp1 (:), dp2(:,:)
    real(ESMF_KIND_R8), pointer :: dpf1(:), dpf2(:,:)
    real(ESMF_KIND_R8), pointer :: dpw1(:), dpw2(:,:)
    integer                     :: nf, n1, n2
    integer                     :: cnt, cntf
    character(ESMF_MAXSTR)      :: FBoutfld, FBoutpre, FBoutname, FBfld, FBname, FBpre
    integer                     :: lrank
    logical                     :: ldocument, match
    character(len=128)          :: lstring
    character(len=*),parameter  :: subname='(med_merge_fbx)'

    !    if (dbug_flag > 5) then
    !      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    !    endif
    rc = ESMF_SUCCESS

    ldocument = .false.
    if (present(document)) ldocument=document

    lstring = ""
    if (present(string)) lstring=string

    call shr_nuopc_methods_FB_getNameN(FBout, n, FBoutfld, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    n2 = scan(FBoutfld,'_')-1
    FBoutpre  = trim(FBoutfld(1:n2))
    FBoutname = trim(FBoutfld(scan(FBoutfld,'_'):))
    call shr_nuopc_methods_FB_GetFldPtr(FBout, trim(FBoutfld), fldptr1=dp1, fldptr2=dp2, rank=lrank, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !    write(tmpstr,*) subname,'tcx FBoutfld=',n,trim(FBoutfld)
    !    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    if (present(FBw) .and. present(fldw)) then
      call shr_nuopc_methods_FB_GetFldPtr(FBw, trim(fldw), fldptr1=dpw1, fldptr2=dpw2, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call ESMF_FieldBundleGet(FB, fieldCount=cntf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    !    write(tmpstr,*) subname,'tcx cntf=',cntf
    !    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    do nf = 1,cntf
      call shr_nuopc_methods_FB_getNameN(FB, nf, FBfld, rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      n2 = scan(FBfld,'_')-1
      FBpre  = trim(FBfld(1:n2))
      FBname = trim(FBfld(scan(FBfld,'_'):))

      if (trim(FBfld) == trim(FBoutfld)) then
        ! this is just a copy of identical full field names

        call shr_nuopc_methods_FB_GetFldPtr(FB, trim(FBfld), dpf1, dpf2, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (document) call ESMF_LogWrite(trim(subname)//":"//trim(lstring)//":"//trim(FBoutfld)//" ="//trim(FBfld), ESMF_LOGMSG_INFO, rc=dbrc)
        if (lrank == 1) then
          dp1(:) = dpf1(:)
        else
          dp2(:,:) = dpf2(:,:)
        endif

      elseif (index(FBoutpre,'x') > 0 .and. trim(FBname) == trim(FBoutname)) then
        ! this is a merge, FBoutpre must have an x in it.

        ! this checks whether a character in FBoutpre matches a character in FBpre that is NOT F or x.
        ! if so, then this term should be merged.  For instance, if FBoutpre=Faxx and FBpre is Faii then the
        ! "a" matches and a merge is fine.  If FBoutpre=Faxx and FBpre is Fioi then there is not match and
        ! that term should not be merged.  This only happens for "F" merges.  If FBpre has an x in it, then
        ! it can always be merged.
        match = .true.
        if (index(FBoutpre,'F') > 0) match = .false.
        if (index(FBpre,'x') > 0) match = .true.
        n1 = 0
        do while (n1 < len_trim(FBoutpre) .and. .not.match)
          n1 = n1 + 1
          if (scan(FBoutpre(n1:n1),'Fx') == 0 .and. scan(FBpre,FBoutpre(n1:n1)) > 0) match = .true.
        enddo

        if (match) then
          call shr_nuopc_methods_FB_GetFldPtr(FB, trim(FBfld), dpf1, dpf2, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          if (present(FBw) .and. present(fldw)) then
            if (document) call ESMF_LogWrite(trim(subname)//":"//trim(lstring)//":"//trim(FBoutfld)//"+="//trim(FBfld)//"*"//trim(fldw), ESMF_LOGMSG_INFO, rc=dbrc)
            if (lrank == 1) then
              dp1(:) = dp1(:) + dpf1(:)*dpw1(:)
            else
              dp2 = dp2 + dpf2*dpw2
            endif
          else
            if (document) call ESMF_LogWrite(trim(subname)//":"//trim(lstring)//":"//trim(FBoutfld)//"+="//trim(FBfld), ESMF_LOGMSG_INFO, rc=dbrc)
            if (lrank == 1) then
              dp1 = dp1 + dpf1
            else
              dp2 = dp2 + dpf2
            endif
          endif
        endif
      endif
    enddo

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    ! if (dbug_flag > 5) then
    !    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    ! endif

  end subroutine med_merge_fbx

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

end module med_merge_mod
