module med_phases_prep_ocn_mod

  !-----------------------------------------------------------------------------
  ! Carry out fast accumulation for the ocean
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use shr_kind_mod            , only : CL=>SHR_KIND_CL, CS=>SHR_KIND_CS
  use esmFlds                 , only : compatm, compocn, compice, compglc, comprof, ncomps, compname 
  use esmFlds                 , only : fldListFr, fldListTo
  use esmFlds                 , only : fldListMed_aoflux_o
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_accum
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_FldChk
  use med_constants_mod       , only : med_constants_dbug_flag
  use med_merge_mod           , only : med_merge_auto
  use med_map_mod             , only : med_map_FB_Regrid_Norm 
  use med_internalstate_mod   , only : InternalState, logunit

  implicit none
  private

  integer           , parameter :: dbug_flag = med_constants_dbug_flag
  character(*)      , parameter :: u_FILE_u  = __FILE__
  integer                       :: dbrc
  logical                       :: mastertask

  ! TODO: the calculation needs to be set at run time based on receiving it from the ocean
  real(ESMF_KIND_R8)            :: flux_epbalfact = 1._ESMF_KIND_R8 

  public :: med_phases_prep_ocn_map
  public :: med_phases_prep_ocn_merge
  public :: med_phases_prep_ocn_accum_fast
  public :: med_phases_prep_ocn_accum_avg

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_map(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)         :: is_local
    integer                     :: n1, ncnt
    character(len=*), parameter :: subname='(med_phases_prep_ocn_map)'
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

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compocn), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ncnt == 0) then
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compocn), returning", &
               ESMF_LOGMSG_INFO, rc=dbrc)
       endif
       RETURN
    end if

    !---------------------------------------
    !--- map all fields in FBImp that have 
    !--- active ocean coupling to the ocean grid
    !---------------------------------------

    do n1 = 1,ncomps
       if (is_local%wrap%med_coupling_active(n1,compocn)) then
          call med_map_FB_Regrid_Norm( &
               fldListFr(n1)%flds, n1, compocn, &
               is_local%wrap%FBImp(n1,n1), &
               is_local%wrap%FBImp(n1,compocn), &
               is_local%wrap%FBFrac(n1), &
               is_local%wrap%FBNormOne(n1,compocn,:), &
               is_local%wrap%RH(n1,compocn,:), &
               string=trim(compname(n1))//'2'//trim(compname(compocn)), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    enddo

  end subroutine med_phases_prep_ocn_map

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_merge(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    logical,save                :: first_call = .true.
    type(InternalState)         :: is_local
    integer                     :: n, n1, ncnt
    integer                     :: lsize
    real(ESMF_KIND_R8), pointer :: dataptr1(:)
    real(ESMF_KIND_R8), pointer :: ifrac(:), ofrac(:)
    real(ESMF_KIND_R8), pointer :: ifracr(:), ofracr(:)
    real(ESMF_KIND_R8), pointer :: avsdr(:), avsdf(:)
    real(ESMF_KIND_R8), pointer :: anidr(:), anidf(:)
    real(ESMF_KIND_R8), pointer :: swvdf(:), swndf(:)
    real(ESMF_KIND_R8), pointer :: swvdr(:), swndr(:)
    real(ESMF_KIND_R8), pointer :: swpen(:), swnet(:)
    real(ESMF_KIND_R8)          :: ifrac_scaled, ofrac_scaled
    real(ESMF_KIND_R8)          :: ifracr_scaled, ofracr_scaled
    real(ESMF_KIND_R8)          :: frac_sum
    real(ESMF_KIND_R8)          :: fswabsv, fswabsi
    character(len=*), parameter :: ice_fraction_name = 'Si_ifrac'
    character(len=*), parameter :: subname='(med_phases_prep_ocn_merge)'
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

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compocn), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt == 0) then
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compocn), returning", &
               ESMF_LOGMSG_INFO, rc=dbrc)
       endif
       RETURN
    end if

    !---------------------------------------
    !--- auto merges to ocn
    !---------------------------------------

    call med_merge_auto(trim(compname(compocn)), &
         is_local%wrap%FBExp(compocn), is_local%wrap%FBFrac(compocn), &
         is_local%wrap%FBImp(:,compocn), fldListTo(compocn), &
         FBMed1=is_local%wrap%FBMed_aoflux_o, &
         document=first_call, string='(merge_to_ocn)', mastertask=mastertask, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------

    if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_rain', rc=rc)) then
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rain' , dataptr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1(:) = dataptr1(:) * flux_epbalfact 
       write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_rain by flux_epbalfact '
    end if
    if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_snow', rc=rc)) then
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_snow' , dataptr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1(:) = dataptr1(:) * flux_epbalfact 
       write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_snow by flux_epbalfact '
    end if
    if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_prec', rc=rc)) then
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_prec' , dataptr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1(:) = dataptr1(:) * flux_epbalfact 
       write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_prec by flux_epbalfact '
    end if
    if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_rofl', rc=rc)) then
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofl' , dataptr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1(:) = dataptr1(:) * flux_epbalfact 
       write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_rofl by flux_epbalfact '
    end if
    if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_rofi', rc=rc)) then
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofi' , dataptr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1(:) = dataptr1(:) * flux_epbalfact 
       write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_rofi by flux_epbalfact '
    end if

    !---------------
    ! Compute swnet to send to ocean
    !---------------
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac' , ifrac, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac' , ofrac, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrad' , ifracr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrad' , ofracr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdr' , avsdr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidr' , anidr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdf' , avsdf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidf' , anidf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdr', swvdr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndr', swndr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdf', swvdf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndf', swndf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (is_local%wrap%comp_present(compice)) then
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen', swpen, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet',  swnet, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,size(ofrac)
       ifrac_scaled = ifrac(n)
       ofrac_scaled = ofrac(n)
       frac_sum = ifrac(n) + ofrac(n)
       if (frac_sum /= 0._ESMF_KIND_r8) then
          ifrac_scaled = ifrac(n) / (frac_sum)
          ofrac_scaled = ofrac(n) / (frac_sum)
       endif
       ifracr_scaled = ifracr(n)
       ofracr_scaled = ofracr(n)
       frac_sum = ifracr(n) + ofracr(n)
       if (frac_sum /= 0._ESMF_KIND_r8) then
          ifracr_scaled = ifracr(n) / (frac_sum)
          ofracr_scaled = ofracr(n) / (frac_sum)
       endif
       fswabsv =  swvdr(n) * (1.0_ESMF_KIND_R8 - avsdr(n)) + &
                  swvdf(n) * (1.0_ESMF_KIND_R8 - avsdf(n))
       fswabsi =  swndr(n) * (1.0_ESMF_KIND_R8 - anidr(n)) + &
                  swndf(n) * (1.0_ESMF_KIND_R8 - anidf(n))
       if (is_local%wrap%comp_present(compice)) then
          swnet(n) = ofracr_scaled*(fswabsv + fswabsi) + ifrac_scaled*swpen(n)
       else
          swnet(n) = ofracr_scaled*(fswabsv + fswabsi)
       end if
       ! if (i2o_per_cat) then
       !   Sf_ofrac(n)  = ofrac(n)
       !   Sf_ofracr(n) = ofracr(n)
       !   Foxx_swnet_ofracr(n) = (fswabsv + fswabsi) * ofracr_scaled
       ! end if
    end do
    ! TODO: document merging

#if (1 == 0)
    ! atm and ice fraction
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), trim(ice_fraction_name), icewgt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    lb1 = lbound(icewgt,1); ub1 = ubound(icewgt,1)
    lb2 = lbound(icewgt,2); ub2 = ubound(icewgt,2)
    allocate(atmwgt(lb1:ub1,lb2:ub2), customwgt(lb1:ub1,lb2:ub2))
    do j = lb2,ub2
       do i = ub1,ib1
          atmwgt(i,j) = 1.0_ESMF_KIND_R8 - icewgt(i,j)
       enddo
    enddo

    !-------------
    ! mean_evap_rate = mean_laten_heat_flux * (1-ice_fraction)/const_lhvap
    !-------------
    !    customwgt = atmwgt / const_lhvap
    !    call shr_nuopc_methods_FB_FieldMerge(is_local%wrap%FBExp(compocn),'mean_evap_rate' , &
    !         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='mean_laten_heat_flux', wgtA=customwgt, rc=rc)
    !    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------
    ! field_for_ocn = field_from_atm * (1-ice_fraction)
    !-------------
    call shr_nuopc_methods_FB_FieldMerge(is_local%wrap%FBExp(compocn),'mean_fprec_rate' , &
         FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='mean_fprec_rate', wgtA=atmwgt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------
    ! End merges
    !-------------
    deallocate(atmwgt,customwgt)

    if (dbug_flag > 1) then
       call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compocn), string=trim(subname)//' FB4ocn_AFmrg ', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
#endif

    if (dbug_flag > 1) then
       call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compocn), string=trim(subname)//' FBexp(compocn) ', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    first_call = .false.

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_prep_ocn_merge

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_accum_fast(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Carry out fast accumulation for the ocean

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(InternalState)         :: is_local
    integer                     :: i,j,n,n1,ncnt
    character(len=*), parameter :: subname='(med_phases_accum_fast)'
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

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compocn), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt == 0) then
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compocn), returning", &
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
    !--- ocean accumulator
    !---------------------------------------

    call shr_nuopc_methods_FB_accum(is_local%wrap%FBExpAccum(compocn), is_local%wrap%FBExp(compocn), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    is_local%wrap%FBExpAccumCnt(compocn) = is_local%wrap%FBExpAccumCnt(compocn) + 1

    if (dbug_flag > 1) then
       call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExpAccum(compocn), &
            string=trim(subname)//' FBaccOcn_AFaccum ', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_prep_ocn_accum_fast

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_accum_avg(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Prepares the OCN import Fields.

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(InternalState)         :: is_local
    integer                     :: i,j,n,n1,ncnt
    character(len=*),parameter  :: subname='(med_phases_prep_ocn_accum_avg)'
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

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compocn), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt == 0) then
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compocn), returning", &
               ESMF_LOGMSG_INFO, rc=dbrc)
       endif
       RETURN
    end if

  end subroutine med_phases_prep_ocn_accum_avg

end module med_phases_prep_ocn_mod
