module med_phases_prep_wav_mod

  !-----------------------------------------------------------------------------
  ! Mediator Phases
  !-----------------------------------------------------------------------------

  implicit none
  private

  character(*)      , parameter :: u_FILE_u  = __FILE__

  public  :: med_phases_prep_wav

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

    subroutine med_phases_prep_wav(gcomp, rc)
      use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
      use ESMF, only : ESMF_GridComp, ESMF_Clock, ESMF_Time
      use ESMF, only : ESMF_GridCompGet, ESMF_FieldBundleGet, ESMF_ClockGet, ESMF_TimeGet
      use ESMF, only : ESMF_ClockPrint
      use med_constants_mod, only : CS
      use esmFlds                 , only : compwav, ncomps, compname
      use esmFlds                 , only : fldListFr, fldListTo
      use shr_nuopc_methods_mod   , only : shr_nuopc_methods_ChkErr
      use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_diagnose
      use med_constants_mod       , only : dbug_flag=>med_constants_dbug_flag
      use med_merge_mod           , only : med_merge_auto
      use med_map_mod             , only : med_map_FB_Regrid_Norm
      use med_internalstate_mod   , only : InternalState, mastertask
      use perf_mod                , only : t_startf, t_stopf
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc

      ! Prepares the WAV import Fields.

      ! local variables
      type(ESMF_Clock)            :: clock
      type(ESMF_Time)             :: time
      character(len=CS)           :: timestr
      type(InternalState)         :: is_local
      integer                     :: i,j,n,n1,ncnt
      logical,save                :: first_call = .true.
      integer :: dbrc
      character(len=*),parameter  :: subname='(med_phases_prep_wav)'
      !---------------------------------------
      call t_startf('MED:'//subname)
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
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

      call ESMF_FieldBundleGet(is_local%wrap%FBExp(compwav), fieldCount=ncnt, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (ncnt == 0) then
         if (dbug_flag > 5) then
            call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compwav), returning", &
                 ESMF_LOGMSG_INFO, rc=dbrc)
         endif
      else

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
         !--- map to create FBimp(:,compwav)
         !---------------------------------------

         do n1 = 1,ncomps
            if (is_local%wrap%med_coupling_active(n1,compwav)) then
               call med_map_FB_Regrid_Norm( &
                    fldListFr(n1)%flds, n1, compwav, &
                    is_local%wrap%FBImp(n1,n1), &
                    is_local%wrap%FBImp(n1,compwav), &
                    is_local%wrap%FBFrac(n1), &
                    is_local%wrap%FBNormOne(n1,compwav,:), &
                    is_local%wrap%RH(n1,compwav,:), &
                    string=trim(compname(n1))//'2'//trim(compname(compwav)), rc=rc)
               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
            endif
         enddo

         !---------------------------------------
         !--- auto merges
         !---------------------------------------

         call med_merge_auto(trim(compname(compwav)), &
              is_local%wrap%FBExp(compwav), is_local%wrap%FBFrac(compwav), &
              is_local%wrap%FBImp(:,compwav), fldListTo(compwav), &
              document=first_call, string='(merge_to_wav)', mastertask=mastertask, rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

         if (dbug_flag > 1) then
            call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compwav), string=trim(subname)//' FBexp(compwav) ', rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         endif

         !---------------------------------------
         !--- custom calculations
         !---------------------------------------

         !---------------------------------------
         !--- update local scalar data
         !---------------------------------------

         !is_local%wrap%scalar_data(1) =

         !---------------------------------------
         !--- clean up
         !---------------------------------------

         first_call = .false.
      endif
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
      call t_stopf('MED:'//subname)

    end subroutine med_phases_prep_wav

end module med_phases_prep_wav_mod
