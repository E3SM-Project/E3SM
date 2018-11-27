module med_phases_prep_ice_mod

  !-----------------------------------------------------------------------------
  ! Mediator Phases
  !-----------------------------------------------------------------------------

  implicit none
  private

  character(*)      , parameter :: u_FILE_u  = __FILE__

  public  :: med_phases_prep_ice

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_ice(gcomp, rc)

    ! Prepares the ICE import Fields.

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_ClockPrint
    use ESMF                  , only : ESMF_FieldBundleGet
    use med_constants_mod     , only : CL, CS, R8
    use esmFlds               , only : compatm, compice, comprof, compglc, ncomps, compname
    use esmFlds               , only : fldListFr, fldListTo
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FldChk
    use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
    use med_merge_mod         , only : med_merge_auto
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use med_internalstate_mod , only : InternalState, logunit, mastertask
    use perf_mod              , only : t_startf, t_stopf
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: clock
    type(ESMF_Time)            :: time
    character(len=64)          :: timestr
    type(InternalState)        :: is_local
    real(R8), pointer          :: dataPtr1(:), dataPtr2(:), dataPtr3(:), dataPtr4(:)
    integer                    :: i,n,n1,ncnt
    character(len=CS)          :: fldname
    real(R8), pointer          :: dataptr(:)
    character(len=1024)        :: msgString
    integer                    :: dbrc
    ! TODO: the calculation needs to be set at run time based on receiving it from the ocean
    real(R8)                   :: flux_epbalfact = 1._R8
    logical,save               :: first_call = .true.
    character(len=*),parameter :: subname='(med_phases_prep_ice)'
    !---------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compice), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt == 0) then
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compice), returning", &
               ESMF_LOGMSG_INFO, rc=dbrc)
       endif
       RETURN
    end if

    !---------------------------------------
    !--- Get the current time from the clock
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, clock=clock)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(time,timestring=timestr)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif
#if DEBUG
    if (mastertask) then
       call ESMF_ClockPrint(clock, options="currTime", preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
#endif
    !---------------------------------------
    !--- map to create FBimp(:,compice)
    !---------------------------------------

    do n1 = 1,ncomps
       if (is_local%wrap%med_coupling_active(n1,compice)) then
          call med_map_FB_Regrid_Norm( &
               fldListFr(n1)%flds, n1, compice, &
               is_local%wrap%FBImp(n1,n1), &
               is_local%wrap%FBImp(n1,compice), &
               is_local%wrap%FBFrac(n1), &
               is_local%wrap%FBNormOne(n1,compice,:), &
               is_local%wrap%RH(n1,compice,:), &
               string=trim(compname(n1))//'2'//trim(compname(compice)), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    enddo

    !---------------------------------------
    !--- auto merges
    !---------------------------------------

    call med_merge_auto(trim(compname(compice)), &
         is_local%wrap%FBExp(compice), is_local%wrap%FBFrac(compice), &
         is_local%wrap%FBImp(:,compice), fldListTo(compice), &
         document=first_call, string='(merge_to_ice)', mastertask=mastertask, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compice), string=trim(subname)//' FBexp(compice) ', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------

    if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compice), 'Faxa_rain', rc=rc)) then
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compice), 'Faxa_rain' , dataptr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1(:) = dataptr1(:) * flux_epbalfact
       if (first_call .and. mastertask) then
          write(logunit,'(a)')'(merge_to_ice): Scaling Faxa_rain by flux_epbalfact '
       end if
    end if
    if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compice), 'Faxa_snow', rc=rc)) then
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compice), 'Faxa_snow' , dataptr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1(:) = dataptr1(:) * flux_epbalfact
       if (first_call .and. mastertask) then
          write(logunit,'(a)')'(merge_to_ice): Scaling Faxa_snow by flux_epbalfact '
       end if
    end if
    if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compice), 'Fixx_rofi', rc=rc)) then
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compice), 'Fixx_rofi' , dataptr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1(:) = dataptr1(:) * flux_epbalfact
       if (first_call .and. mastertask) then
          write(logunit,'(a)')'(merge_to_ice): Scaling Fixx_rofi by flux_epbalfact '
       end if
    end if

    !---------------------------------------
    !--- update local scalar data
    !---------------------------------------

    !is_local%wrap%scalar_data(1) =

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ice

end module med_phases_prep_ice_mod
