module med_phases_prep_ice_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing ice export from mediator
  !-----------------------------------------------------------------------------

  implicit none
  private

  character(*)      , parameter :: u_FILE_u  = __FILE__

  public  :: med_phases_prep_ice

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_ice(gcomp, rc)

    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet 
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FieldBundleGet, ESMF_RouteHandleIsCreated
    use ESMF                  , only : ESMF_LOGMSG_ERROR, ESMF_FAILURE
    use NUOPC                 , only : NUOPC_IsConnected
    use esmFlds               , only : compatm, compice, comprof, compglc, ncomps, compname
    use esmFlds               , only : fldListFr, fldListTo
    use esmFlds               , only : mapbilnr
    use shr_nuopc_methods_mod , only : fldchk => shr_nuopc_methods_FB_FldChk
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : FB_GetFldPtr => shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FieldRegrid
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getNumFlds
    use med_constants_mod     , only : CS, R8
    use med_infodata_mod      , only : med_infodata, med_infodata_GetData
    use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
    use med_merge_mod         , only : med_merge_auto
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use med_internalstate_mod , only : InternalState, logunit, mastertask
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=64)              :: timestr
    type(InternalState)            :: is_local
    integer                        :: i,n,n1,ncnt
    character(len=CS)              :: fldname
    real(R8), pointer              :: dataptr(:)
    real(R8), pointer              :: temperature(:)
    real(R8), pointer              :: pressure(:)
    real(R8), pointer              :: humidity(:)
    real(R8), pointer              :: air_density(:)
    real(R8), pointer              :: pot_temp(:)
    real(R8)                       :: precip_fact
    character(len=CS)              :: cvalue
    integer                        :: dbrc
    character(len=64), allocatable :: fldnames(:)
    logical                        :: first_precip_fact_call = .true.
    character(len=*),parameter     :: subname='(med_phases_prep_ice)'
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

    call shr_nuopc_methods_FB_getNumFlds(is_local%wrap%FBExp(compice), trim(subname)//"FBexp(compice)", ncnt, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       !--- map to create FBImp(:,compice)
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
       !--- auto merges to create FBExp(compice)
       !---------------------------------------

       call med_merge_auto(trim(compname(compice)), &
            is_local%wrap%FBExp(compice), is_local%wrap%FBFrac(compice), &
            is_local%wrap%FBImp(:,compice), fldListTo(compice), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       !--- custom calculations
       !---------------------------------------

       ! application of precipitation factor from ocean

       ! TODO (mvertens, 2019-03-18): precip_fact here is not valid if
       ! the component does not send it - so need to check if it is a
       ! valid value (determine the ranges)
       ! For now - hardwire it to 1 until this is resolved
       ! call med_infodata_GetData(med_infodata, precip_fact=precip_fact, rc=rc)
       ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       precip_fact = 1.0_R8

       if (precip_fact /= 1.0_R8) then
          if (first_precip_fact_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ice): Scaling rain, snow, liquid and ice runoff by precip_fact '
             first_precip_fact_call = .false.
          end if
          write(cvalue,*) precip_fact
          call ESMF_LogWrite(trim(subname)//" precip_fact is "//trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

          allocate(fldnames(3))
          fldnames = (/'Faxa_rain', 'Faxa_snow', 'Fixx_rofi'/)
          do n = 1,size(fldnames)
             if (fldchk(is_local%wrap%FBExp(compice), trim(fldnames(n)), rc=rc)) then
                call FB_GetFldPtr(is_local%wrap%FBExp(compice), trim(fldnames(n)) , dataptr, rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                dataptr(:) = dataptr(:) * precip_fact
             end if
          end do
          deallocate(fldnames)
       end if

       ! If either air density or ptem from atm is not available - then need pbot since it will be
       ! required for either calculation
       if ( .not. fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_dens',rc=rc) .or. &
            .not. fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem',rc=rc)) then 

          ! Determine Sa_pbot on the ice grid and get a pointer to it
          if (.not. fldchk(is_local%wrap%FBExp(compice), 'Sa_pbot',rc=rc)) then
             if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,compice,mapbilnr))) then
                call ESMF_LogWrite(trim(subname)//": ERROR bilinr RH not available for atm->ice", &
                     ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
                rc = ESMF_FAILURE
                return
             end if
             call shr_nuopc_methods_FB_FieldRegrid( &
                  is_local%wrap%FBImp(compatm,compatm), 'Sa_pbot', &
                  is_local%wrap%FBImp(compatm,compice), 'Sa_pbot', &
                  is_local%wrap%RH(compatm,compice,mapbilnr), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compice), 'Sa_pbot', pressure, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Get a pointer to Sa_tbot on the ice grid
          call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compice), 'Sa_tbot', temperature, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! compute air density as a custom calculation
       if ( .not. fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_dens',rc=rc)) then
          call ESMF_LogWrite(trim(subname)//": computing air density as a custom calculation", ESMF_LOGMSG_INFO, rc=dbrc)

          call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compice), 'Sa_shum', humidity, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBExp(compice), 'Sa_dens', air_density, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          do n = 1,size(temperature)
             if (temperature(n) /= 0._R8) then
                air_density(n) = pressure(n) / (287.058_R8*(1._R8 + 0.608_R8*humidity(n))*temperature(n))
             else
                air_density(n) = 0._R8
             endif
          end do
       end if

       ! compute potential temperature as a custom calculation
       if (.not. fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem',rc=rc)) then
          call ESMF_LogWrite(trim(subname)//": computing potential temp as a custom calculation", ESMF_LOGMSG_INFO, rc=dbrc)

          call FB_GetFldPtr(is_local%wrap%FBExp(compice), 'Sa_ptem', pot_temp, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          do n = 1,size(temperature)
             if (pressure(n) /= 0._R8) then
                pot_temp(n) = temperature(n) * (100000._R8/pressure(n))**0.286_R8 ! Potential temperature (K)
             else
                pot_temp(n) = 0._R8
             end if
          end do
       end if

       !if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compice), string=trim(subname)//' FBexp(compice) ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       !end if

       !---------------------------------------
       !--- update local scalar data
       !---------------------------------------

       !is_local%wrap%scalar_data(1) =

       !---------------------------------------
       !--- clean up
       !---------------------------------------

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ice

end module med_phases_prep_ice_mod
