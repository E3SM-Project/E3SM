module med_phases_prep_lnd_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing land export from mediator
  !-----------------------------------------------------------------------------

  implicit none
  private

  character(*) , parameter :: u_FILE_u = &
       __FILE__

  public  :: med_phases_prep_lnd

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_lnd(gcomp, rc)

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_ClockPrint
    use ESMF                  , only : ESMF_FieldBundleGet
    use esmFlds               , only : complnd, compatm, ncomps, compname 
    use esmFlds               , only : fldListFr, fldListTo
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getNumFlds
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
    use shr_nuopc_scalars_mod , only : flds_scalar_name, flds_scalar_num
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nextsw_cday
    use med_constants_mod     , only : R8, dbug_flag=>med_constants_dbug_flag
    use med_merge_mod         , only : med_merge_auto
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use med_internalstate_mod , only : InternalState
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n1,ncnt
    real(r8)            :: nextsw_cday
    character(len=*), parameter :: subname='(med_phases_prep_lnd)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

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

    call shr_nuopc_methods_FB_getNumFlds(is_local%wrap%FBExp(complnd), trim(subname)//"FBexp(complnd)", ncnt, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       !--- map to create FBimp(:,complnd)
       !---------------------------------------

       do n1 = 1,ncomps
          if (is_local%wrap%med_coupling_active(n1,complnd)) then
             call med_map_FB_Regrid_Norm( &
                  fldListFr(n1)%flds, n1, complnd, &
                  is_local%wrap%FBImp(n1,n1), &
                  is_local%wrap%FBImp(n1,complnd), &
                  is_local%wrap%FBFrac(n1), &
                  is_local%wrap%FBFrac(complnd), &
                  is_local%wrap%FBNormOne(n1,complnd,:), &
                  is_local%wrap%RH(n1,complnd,:), &
                  string=trim(compname(n1))//'2'//trim(compname(complnd)), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       enddo

       !---------------------------------------
       !--- auto merges to create FBExp(complnd)
       !---------------------------------------

       call med_merge_auto(trim(compname(complnd)), &
            is_local%wrap%FBExp(complnd), &
            is_local%wrap%FBFrac(complnd), &
            is_local%wrap%FBImp(:,complnd), &
            fldListTo(complnd), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(complnd), &
               string=trim(subname)//' FBexp(complnd) ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- custom calculations
       !---------------------------------------

       !---------------------------------------
       !--- update scalar data
       !---------------------------------------

       ! send nextsw_cday to land - first obtain it from atm import 
       call shr_nuopc_methods_State_GetScalar(&
            scalar_value=nextsw_cday, scalar_id=flds_scalar_index_nextsw_cday, &
            state=is_local%wrap%NstateImp(compatm), flds_scalar_name=flds_scalar_name, &
            flds_scalar_num=flds_scalar_num, rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_State_SetScalar(&
            scalar_value=nextsw_cday, scalar_id=flds_scalar_index_nextsw_cday, &
            state=is_local%wrap%NstateExp(complnd), flds_scalar_name=flds_scalar_name, &
            flds_scalar_num=flds_scalar_num, rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       !--- clean up
       !---------------------------------------

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_lnd

end module med_phases_prep_lnd_mod
