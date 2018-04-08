module med_phases_accum_fast_mod

  !-----------------------------------------------------------------------------
  ! Mediator Phases
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use shr_kind_mod            , only : CL=>SHR_KIND_CL, CS=>SHR_KIND_CS
  use esmFlds                 , only : compatm, complnd, compocn, compice, comprof, compglc
  use esmFlds                 , only : ncomps, compname 
  use esmFlds                 , only : fldListFr, fldListTo
  use esmFlds                 , only : fldListMed_aoflux_a, fldListMed_aoflux_o
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_reset
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_accum
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

  public  :: med_phases_accum_fast

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

  subroutine med_phases_accum_fast(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Carry out fast accumulation for the ocean

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(InternalState)         :: is_local
    integer                     :: i,j,n,n1,ncnt
    logical,save                :: first_call = .true.
    character(len=*), parameter :: ice_fraction_name = 'Si_ifrac'
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
    !--- map all fields in FBImp that have 
    !--- active ocean coupling to the ocean grid
    !---------------------------------------

    do n1 = 1,ncomps
       if (is_local%wrap%med_coupling_active(n1,compocn)) then
          call med_map_FB_Regrid_Norm( &
               fldListFr(n1)%flds, n1, compocn, &
               is_local%wrap%FBImp(n1,n1), &
               is_local%wrap%FBImp(n1,compocn), &
               is_local%wrap%FBFrac(compocn), &
               is_local%wrap%FBNormOne(n1,compocn,:), &
               is_local%wrap%RH(n1,compocn,:), &
               string=trim(compname(n1))//'2'//trim(compname(compocn)), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    enddo

    !---------------------------------------
    !--- auto merges to ocn
    !---------------------------------------

    call shr_nuopc_methods_FB_reset(is_local%wrap%FBExp(compocn), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call med_merge_auto(is_local%wrap%FBExp(compocn), &
         FB1=is_local%wrap%FBImp(compatm,compocn)   , FB1w=is_local%wrap%FBfrac(compocn), fldw1='afrac', &
         FB2=is_local%wrap%FBMed_aoflux_o           , FB2w=is_local%wrap%FBfrac(compocn), fldw2='afrac', &
         FB3=is_local%wrap%FBImp(compice,compocn)   , FB3w=is_local%wrap%FBfrac(compocn), fldw3='ifrac', &
         FB4=is_local%wrap%FBImp(comprof,compocn)   , &
         FB5=is_local%wrap%FBImp(compglc,compocn)   , &
         document=first_call, string=subname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations to ocn
    !---------------------------------------

#if (1 == 0)
    ! atm and ice fraction
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), trim(ice_fraction_name), icewgt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(atmwgt   (lbound(icewgt,1):ubound(icewgt,1),lbound(icewgt,2):ubound(icewgt,2)))
    allocate(customwgt(lbound(icewgt,1):ubound(icewgt,1),lbound(icewgt,2):ubound(icewgt,2)))
    do j=lbound(icewgt,2),ubound(icewgt,2)
       do i=lbound(icewgt,1),ubound(icewgt,1)
          atmwgt = 1.0_ESMF_KIND_R8 - icewgt
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
       call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compocn), &
            string=trim(subname)//' FB4ocn_AFmrg ', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
#endif

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

    first_call = .false.

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_accum_fast

end module med_phases_accum_fast_mod
