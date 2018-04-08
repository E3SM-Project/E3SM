module med_phases_prep_ice_mod

  !-----------------------------------------------------------------------------
  ! Mediator Phases
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use shr_kind_mod            , only : CL=>SHR_KIND_CL, CS=>SHR_KIND_CS
  use esmFlds                 , only : compatm, compocn, compice, comprof, compglc 
  use esmFlds                 , only : ncomps, compname 
  use esmFlds                 , only : flds_scalar_name
  use esmFlds                 , only : fldListFr, fldListTo
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_reset
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_GetFldPtr
  use med_constants_mod       , only : med_constants_dbug_flag
  use med_constants_mod       , only : med_constants_czero
  use med_merge_mod           , only : med_merge_auto
  use med_map_mod             , only : med_map_FB_Regrid_Norm 
  use med_internalstate_mod   , only : InternalState

  implicit none
  private

  integer           , parameter :: dbug_flag = med_constants_dbug_flag
  real(ESMF_KIND_R8), parameter :: czero     = med_constants_czero
  character(*)      , parameter :: u_FILE_u  = __FILE__
  integer                       :: dbrc
  logical                       :: mastertask

  public  :: med_phases_prep_ice

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_ice(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Prepares the ICE import Fields.

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(InternalState)         :: is_local
    real(ESMF_KIND_R8), pointer :: dataPtr1(:), dataPtr2(:), dataPtr3(:), dataPtr4(:)
    integer                     :: i,n,n1,ncnt
    logical,save                :: first_call = .true.
    character(len=CS)           :: fldname         
    real(ESMF_KIND_R8), pointer :: dataptr(:)      
    character(len=1024)         :: msgString       
    character(len=*),parameter  :: subname='(med_phases_prep_ice)'
    !---------------------------------------

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

    if (mastertask) then
       call ESMF_ClockPrint(clock, options="currTime", preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

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

    call shr_nuopc_methods_FB_reset(is_local%wrap%FBExp(compice), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call med_merge_auto(is_local%wrap%FBExp(compice), &
         FB1=is_local%wrap%FBImp(compocn,compice), FB1w=is_local%wrap%FBfrac(compice), fldw1='ofrac', &
         FB2=is_local%wrap%FBImp(compatm,compice), FB2w=is_local%wrap%FBfrac(compice), fldw2='afrac', &
         FB3=is_local%wrap%FBImp(compglc,compice), &
         FB4=is_local%wrap%FBImp(comprof,compice), &
         document=first_call, string=subname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------

    ! TODO: document custom merges below
    ! TODO: need to obtain flux_epbalfact

    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compice), 'Faxa_rainc', dataPtr1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compice), 'Faxa_rainl', dataPtr2, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compice), 'Faxa_rain' , dataPtr3, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    dataPtr3(:) = dataPtr1(:) + dataPtr2(:)
#if (1 == 0)
    dataPtr3(:) = dataPtr3(:) * flux_epbalfact
#endif

    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compice), 'Faxa_snowc', dataPtr1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compice), 'Faxa_snowl', dataPtr2, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compice), 'Faxa_snow' , dataPtr3, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    dataPtr3(:) = dataPtr1(:) + dataPtr2(:)
#if (1 == 0)
    dataPtr3(:) = dataPtr3(:) * flux_epbalfact
#endif

#if (1 == 0)
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compglc,compice), 'Figg_rofi', dataPtr1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(comprof,compice), 'Firr_rofi', dataPtr2, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compice), 'Fixx_rofi', dataPtr3, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    dataPtr3(:) = dataPtr1(:) + dataPtr2(:)
    dataPtr3(:) = dataPtr3(:) * flux_epbalfact
#endif

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

  end subroutine med_phases_prep_ice

end module med_phases_prep_ice_mod
