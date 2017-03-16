!-----------------------------------------------------------------------
!
MODULE MOSART_physics_mod
! Description: core code of MOSART. Can be incoporated within any land model via a interface module
! 
! Developed by Hongyi Li, 12/29/2011. 
! REVISION HISTORY:
! Jan 2012, only consider land surface water routing, no parallel computation
! May 2012, modified to be coupled with CLM
!-----------------------------------------------------------------------

! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use shr_sys_mod   , only : shr_sys_abort
  use RtmVar        , only : iulog, barrier_timers, wrmflag
  use RunoffMod     , only : Tctl, TUnit, TRunoff, TPara, rtmCTL, &
                             SMatP_eroutUp, avsrc_eroutUp, avdst_eroutUp
  use RtmSpmd       , only : masterproc, mpicom_rof, iam
  use RtmTimeManager, only : get_curr_date, is_new_month
#ifdef INCLUDE_WRM
  use WRM_type_mod  , only : ctlSubwWRM, WRMUnit, StorWater
  use WRM_modules   , only : irrigationExtractionSubNetwork, &
                             irrigationExtractionMainChannel, &
                             Regulation, ExtractionRegulatedFlow
  use WRM_returnflow, only : insert_returnflow_channel, &
                             insert_returnflow_soilcolumn, &
                             estimate_returnflow_deficit
  use WRM_subw_io_mod, only : WRM_readDemand, WRM_computeRelease
#endif
  use rof_cpl_indices, only : nt_rtm, rtm_tracers, nt_nliq
  use perf_mod, only: t_startf, t_stopf
  use mct_mod

  implicit none
  private

  real(r8), parameter :: TINYVALUE = 1.0e-14_r8  ! double precision variable has a significance of about 16 decimal digits
    integer  :: nt               ! loop indices
  real(r8), parameter :: SLOPE1def = 0.1_r8        ! here give it a small value in order to avoid the abrupt change of hydraulic radidus etc.
  real(r8) :: sinatanSLOPE1defr   ! 1.0/sin(atan(slope1))
  public Euler
  public updatestate_hillslope
  public updatestate_subnetwork
  public updatestate_mainchannel

!-----------------------------------------------------------------------
                   
! !PUBLIC MEMBER FUNCTIONS:
  contains

!-----------------------------------------------------------------------
  subroutine Euler
  ! !DESCRIPTION: solve the ODEs with Euler algorithm
    implicit none    
    
    integer :: iunit, idam, m, k, unitUp, cnt, ier   !local index
    real(r8) :: temp_erout, localDeltaT
    real(r8) :: negchan
    integer  :: yr,mon,day,tod
    character(len=*),parameter :: subname = '(Euler)'
    !------------------

    call get_curr_date(yr, mon, day, tod)

    !------------------
    ! WRM prep
    !------------------

#ifdef INCLUDE_WRM
    if (wrmflag) then
       if ( ctlSubwWRM%ReturnFlowFlag > 0) then
          call insert_returnflow_soilcolumn
       endif

       !call readPotentialEvap(trim(theTime))
       if ( is_new_month() ) then
          if (masterproc) write(iulog,*) trim(subname),' updating monthly data at ',yr,mon,day,tod

          ! presently the demand is hardcoded at a monthly time step else need to
          ! be moved out of this if loopon first time step of the month
          call WRM_readDemand()
          call WRM_computeRelease()

       end if
       StorWater%demand = StorWater%demand0 * Tctl%DeltaT
       !supply is set to zero in RtmMod so it can be accumulated there for the budget
       !StorWater%supply = 0._r8
       StorWater%deficit =0._r8
    endif
#endif

    !------------------
    ! hillslope
    !------------------

    call t_startf('mosartr_hillslope')
    do nt=1,nt_rtm
    if (TUnit%euler_calc(nt)) then
    do iunit=rtmCTL%begr,rtmCTL%endr
       if(TUnit%mask(iunit) > 0) then
          call hillslopeRouting(iunit,nt,Tctl%DeltaT)
          TRunoff%wh(iunit,nt) = TRunoff%wh(iunit,nt) + TRunoff%dwh(iunit,nt) * Tctl%DeltaT
          call UpdateState_hillslope(iunit,nt)
! NV WARNNG WARNING JUST FOR TESTING
!          TRunoff%qsub(iunit,nt) = TRunoff%qsub(iunit,nt) / 1.6_r8
          TRunoff%etin(iunit,nt) = (-TRunoff%ehout(iunit,nt) + TRunoff%qsub(iunit,nt)) * TUnit%area(iunit) * TUnit%frac(iunit)
       endif
    end do
    endif
    end do
    call t_stopf('mosartr_hillslope')

!moved inside the subnetwork channel routing and work on wt instead of etin
!#ifdef INCLUDE_WRM
!    if (wrmflag) then
!       call t_startf('mosartr_wrm_IESubN')
!       ! extraction from available surface runoff 
!       if (ctlSubwWRM%ExtractionFlag > 0) then
!          do iunit=rtmCTL%begr,rtmCTL%endr
!             if (TUnit%mask(iunit) > 0) then
!                call irrigationExtractionSubNetwork(iunit, Tctl%DeltaT )
!             endif
!          enddo
!       end if
!       call t_stopf('mosartr_wrm_IESubN')
!    endif
!#endif

    TRunoff%flow = 0._r8
    TRunoff%eroup_lagi = 0._r8
    TRunoff%eroup_lagf = 0._r8
    TRunoff%eroutup_avg = 0._r8
    TRunoff%erlat_avg = 0._r8
    negchan = 9999.0_r8
    do m=1,Tctl%DLevelH2R

       !------------------
       ! subnetwork
       !------------------

       call t_startf('mosartr_subnetwork')    
       TRunoff%erlateral(:,:) = 0._r8
       do nt=1,nt_rtm
       if (TUnit%euler_calc(nt)) then
       do iunit=rtmCTL%begr,rtmCTL%endr
          if(TUnit%mask(iunit) > 0) then
!extraction from subnetwork here from wt
#ifdef INCLUDE_WRM
             if (wrmflag) then
                if (nt == nt_nliq) then
                   if  (ctlSubwWRM%ExtractionFlag > 0 ) then
                      localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
                      call irrigationExtractionSubNetwork(iunit, localDeltaT )
                      call UpdateState_subnetwork(iunit,nt)
                   endif
                endif
             endif
#endif
             localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_t(iunit)
             do k=1,TUnit%numDT_t(iunit)
                call subnetworkRouting(iunit,nt,localDeltaT)
                TRunoff%wt(iunit,nt) = TRunoff%wt(iunit,nt) + TRunoff%dwt(iunit,nt) * localDeltaT
                call UpdateState_subnetwork(iunit,nt)
                TRunoff%erlateral(iunit,nt) = TRunoff%erlateral(iunit,nt)-TRunoff%etout(iunit,nt)
             end do ! numDT_t
             TRunoff%erlateral(iunit,nt) = TRunoff%erlateral(iunit,nt) / TUnit%numDT_t(iunit)
          endif
       end do ! iunit
       endif  ! euler_calc
       end do ! nt
       call t_stopf('mosartr_subnetwork')    

       !------------------
       ! upstream interactions
       !------------------

       if (barrier_timers) then
          call t_startf('mosartr_SMeroutUp_barrier')    
          call mpi_barrier(mpicom_rof,ier)
          call t_stopf('mosartr_SMeroutUp_barrier')    
       endif

       !--- accumulate/average erout at prior timestep (used in eroutUp calc) for budget analysis
       Trunoff%eroup_lagi = Trunoff%eroup_lagi - Trunoff%erout

       call t_startf('mosartr_SMeroutUp')    
       TRunoff%eroutUp = 0._r8
#ifdef NO_MCT
       do iunit=rtmCTL%begr,rtmCTL%endr
       do k=1,TUnit%nUp(iunit)
          unitUp = Tunit%iUp(iunit,k)
          do nt=1,nt_rtm
             TRunoff%eroutUp(iunit,nt) = TRunoff%eroutUp(iunit,nt) + TRunoff%erout(unitUp,nt)
          end do
       end do
       end do
#else
       !--- copy erout into avsrc_eroutUp ---
       call mct_avect_zero(avsrc_eroutUp)
       cnt = 0
       do iunit = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          do nt = 1,nt_rtm
             avsrc_eroutUp%rAttr(nt,cnt) = TRunoff%erout(iunit,nt)
          enddo
       enddo
       call mct_avect_zero(avdst_eroutUp)

       call mct_sMat_avMult(avsrc_eroutUp, sMatP_eroutUp, avdst_eroutUp)

       !--- add mapped eroutUp to TRunoff ---
       cnt = 0
       do iunit = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          do nt = 1,nt_rtm
             TRunoff%eroutUp(iunit,nt) = avdst_eroutUp%rAttr(nt,cnt)
          enddo
       enddo
#endif
       call t_stopf('mosartr_SMeroutUp')    

       TRunoff%eroutup_avg = TRunoff%eroutup_avg + TRunoff%eroutUp
       TRunoff%erlat_avg   = TRunoff%erlat_avg   + TRunoff%erlateral

       !------------------
       ! channel routing
       !------------------

       call t_startf('mosartr_chanroute')    
       do nt=1,nt_rtm
       if (TUnit%euler_calc(nt)) then
       do iunit=rtmCTL%begr,rtmCTL%endr
          if(TUnit%mask(iunit) > 0) then
             localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_r(iunit)
             temp_erout = 0._r8
             do k=1,TUnit%numDT_r(iunit)
                call mainchannelRouting(iunit,nt,localDeltaT)    
                TRunoff%wr(iunit,nt) = TRunoff%wr(iunit,nt) + TRunoff%dwr(iunit,nt) * localDeltaT
! check for negative channel storage
!                if(TRunoff%wr(iunit,1) < -1.e-10) then
!                   write(iulog,*) 'Negative channel storage! ', iunit, TRunoff%wr(iunit,1)
!                   call shr_sys_abort('mosart: negative channel storage')
!                end if
                call UpdateState_mainchannel(iunit,nt)
                temp_erout = temp_erout + TRunoff%erout(iunit,nt) ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
             end do
             temp_erout = temp_erout / TUnit%numDT_r(iunit)
             TRunoff%erout(iunit,nt) = temp_erout
#ifdef INCLUDE_WRM
             if (wrmflag) then
                if (nt == nt_nliq) then
                   localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
                   if (ctlSubwWRM%ExtractionMainChannelFlag > 0 .AND. ctlSubwWRM%ExtractionFlag > 0 ) then
                      call IrrigationExtractionMainChannel(iunit, localDeltaT )
                      if (ctlSubwWRM%TotalDemandFlag > 0 .AND. ctlSubwWRM%ReturnFlowFlag > 0 ) then
                         call insert_returnflow_channel(iunit, localDeltaT )
                      endif
                      ! update main channel storage as well
                      temp_erout = temp_erout - TRunoff%erout(iunit,nt) ! change in erout after regulation and extraction
                      TRunoff%dwr(iunit,nt) =  temp_erout
                      TRunoff%wr(iunit,nt) = TRunoff%wr(iunit,nt) + TRunoff%dwr(iunit,nt) * localDeltaT
                      call UpdateState_mainchannel(iunit,nt)
                   endif
! tcraig, moved out of loop
!                   if ( ctlSubwWRM%RegulationFlag>0 .and. WRMUnit%INVicell(iunit) > 0 .and. WRMUnit%MeanMthFlow(iunit,13) > 0.01_r8 ) then
!                      call Regulation(iunit, localDeltaT)
!                      if ( ctlSubwWRM%ExtractionFlag > 0 ) then
!                         call ExtractionRegulatedFlow(iunit, localDeltaT)
!                      endif
!                   endif
                endif
!                ! do not update wr after regulation or extraction from reservoir release. Because of the regulation, 
!                ! the wr might get to crazy uncontrolled values, assume in this case wr is not changed. The storage in reservoir handles it.
             endif
#endif
             Trunoff%eroup_lagf(iunit,nt) = Trunoff%eroup_lagf(iunit,nt) - Trunoff%erout(iunit,nt)
             TRunoff%flow(iunit,nt) = TRunoff%flow(iunit,nt) - TRunoff%erout(iunit,nt)
          endif
       end do ! iunit
       endif  ! euler_calc
       end do ! nt
       negchan = min(negchan, minval(TRunoff%wr(:,:)))

       call t_stopf('mosartr_chanroute')    
    end do  ! DLevelH2R

! check for negative channel storage
    if (negchan < -1.e-10) then
       write(iulog,*) 'Warning: Negative channel storage found! ',negchan
!       call shr_sys_abort('mosart: negative channel storage')
    endif
    TRunoff%flow = TRunoff%flow / Tctl%DLevelH2R
    TRunoff%eroup_lagi = TRunoff%eroup_lagi / Tctl%DLevelH2R
    TRunoff%eroup_lagf = TRunoff%eroup_lagf / Tctl%DLevelH2R
    TRunoff%eroutup_avg = TRunoff%eroutup_avg / Tctl%DLevelH2R
    TRunoff%erlat_avg = TRunoff%erlat_avg / Tctl%DLevelH2R

    !------------------
    ! WRM Regulation
    ! WRM ExtractionRegulatedFlow
    ! Do not update wr after regulation or extraction from reservoir release. 
    ! Because of the regulation, the wr might get to crazy uncontrolled values, 
    ! assume in this case wr is not changed. The storage in reservoir handles it.
    !------------------

#ifdef INCLUDE_WRM
    if (wrmflag) then
       if (ctlSubwWRM%RegulationFlag>0) then
          ! compute the erowm_reg terms and adjust the flow diagnostic
          do iunit=rtmCTL%begr,rtmCTL%endr
             TRunoff%erowm_regi(iunit,nt_nliq) = -TRunoff%erout(iunit,nt_nliq)
             TRunoff%flow(iunit,nt_nliq) = TRunoff%flow(iunit,nt_nliq) + TRunoff%erout(iunit,nt_nliq)
          enddo
          localDeltaT = Tctl%DeltaT
          call t_startf('mosartr_wrm_Reg')
          do iunit=rtmCTL%begr,rtmCTL%endr
             if (TUnit%mask(iunit) > 0) then
                call Regulation(iunit, localDeltaT)
             endif
          enddo
          call t_stopf('mosartr_wrm_Reg')
          if (ctlSubwWRM%ExtractionFlag > 0 ) then
             call t_startf('mosartr_wrm_ERFlow')
             call ExtractionRegulatedFlow(localDeltaT)
             call t_stopf('mosartr_wrm_ERFlow')
          endif
          !--- now subtract updated erout to update flow calc
          ! compute the erowm_reg terms and adjust the flow diagnostic
          do iunit=rtmCTL%begr,rtmCTL%endr
             TRunoff%erowm_regf(iunit,nt_nliq) = -TRunoff%erout(iunit,nt_nliq)
             TRunoff%flow(iunit,nt_nliq) = TRunoff%flow(iunit,nt_nliq) - TRunoff%erout(iunit,nt_nliq)
          enddo
       endif
    endif
#endif

    !------------------
    ! WRM post Euler updates
    !------------------

#ifdef INCLUDE_WRM
    if (wrmflag) then
       call t_startf('mosartr_wrm_estrfdef')
       call estimate_returnflow_deficit()
       if (ctlSubwWRM%ExtractionFlag > 0) then
          StorWater%deficit = StorWater%demand
       endif
       call t_stopf('mosartr_wrm_estrfdef')
    endif
#endif

  end subroutine Euler

!-----------------------------------------------------------------------

  subroutine hillslopeRouting(iunit, nt, theDeltaT)
  ! !DESCRIPTION: Hillslope routing considering uniform runoff generation across hillslope
    implicit none
    
    integer, intent(in) :: iunit, nt
    real(r8), intent(in) :: theDeltaT    
    character(len=*),parameter :: subname = '(hillslopeRouting)'

!  !TRunoff%ehout(iunit,nt) = -CREHT(TUnit%hslp(iunit), TUnit%nh(iunit), TUnit%Gxr(iunit), TRunoff%yh(iunit,nt))
    TRunoff%ehout(iunit,nt) = -CREHT_nosqrt(TUnit%hslpsqrt(iunit), TUnit%nh(iunit), TUnit%Gxr(iunit), TRunoff%yh(iunit,nt))
    if(TRunoff%ehout(iunit,nt) < 0._r8 .and. &
       TRunoff%wh(iunit,nt) + (TRunoff%qsur(iunit,nt) + TRunoff%ehout(iunit,nt)) * theDeltaT < TINYVALUE) then
         TRunoff%ehout(iunit,nt) = -(TRunoff%qsur(iunit,nt) + TRunoff%wh(iunit,nt) / theDeltaT)
    end if
    TRunoff%dwh(iunit,nt) = (TRunoff%qsur(iunit,nt) + TRunoff%ehout(iunit,nt)) 

  end subroutine hillslopeRouting

!-----------------------------------------------------------------------

  subroutine subnetworkRouting(iunit,nt,theDeltaT)
  ! !DESCRIPTION: subnetwork channel routing
    implicit none    
    integer, intent(in) :: iunit,nt
    real(r8), intent(in) :: theDeltaT
    character(len=*),parameter :: subname = '(subnetworkRouting)'

!  !if(TUnit%tlen(iunit) <= 1e100_r8) then ! if no tributaries, not subnetwork channel routing
    if(TUnit%tlen(iunit) <= TUnit%hlen(iunit)) then ! if no tributaries, not subnetwork channel routing
       TRunoff%etout(iunit,nt) = -TRunoff%etin(iunit,nt)
    else
!     !TRunoff%vt(iunit,nt) = CRVRMAN(TUnit%tslp(iunit), TUnit%nt(iunit), TRunoff%rt(iunit,nt))
       TRunoff%vt(iunit,nt) = CRVRMAN_nosqrt(TUnit%tslpsqrt(iunit), TUnit%nt(iunit), TRunoff%rt(iunit,nt))
       TRunoff%etout(iunit,nt) = -TRunoff%vt(iunit,nt) * TRunoff%mt(iunit,nt)
       if(TRunoff%wt(iunit,nt) + (TRunoff%etin(iunit,nt) + TRunoff%etout(iunit,nt)) * theDeltaT < TINYVALUE) then
          TRunoff%etout(iunit,nt) = -(TRunoff%etin(iunit,nt) + TRunoff%wt(iunit,nt)/theDeltaT)
          if(TRunoff%mt(iunit,nt) > 0._r8) then
             TRunoff%vt(iunit,nt) = -TRunoff%etout(iunit,nt)/TRunoff%mt(iunit,nt)
          end if
       end if
    end if
    TRunoff%dwt(iunit,nt) = TRunoff%etin(iunit,nt) + TRunoff%etout(iunit,nt)

! check stability
!    if(TRunoff%vt(iunit,nt) < -TINYVALUE .or. TRunoff%vt(iunit,nt) > 30) then
!       write(iulog,*) "Numerical error in subnetworkRouting, ", iunit,nt,TRunoff%vt(iunit,nt)
!    end if

  end subroutine subnetworkRouting

!-----------------------------------------------------------------------

  subroutine mainchannelRouting(iunit, nt, theDeltaT)
  ! !DESCRIPTION: main channel routing
    implicit none    
    integer, intent(in) :: iunit, nt
    real(r8), intent(in) :: theDeltaT    
    character(len=*),parameter :: subname = '(mainchannelRouting)'

    if(Tctl%RoutingMethod == 1) then
       call Routing_KW(iunit, nt, theDeltaT)
    else if(Tctl%RoutingMethod == 2) then
       call Routing_MC(iunit, nt, theDeltaT)
    else if(Tctl%RoutingMethod == 3) then
       call Routing_THREW(iunit, nt, theDeltaT)
    else if(Tctl%RoutingMethod == 4) then
       call Routing_DW(iunit, nt, theDeltaT)
    else
       print*, "Please check the routing method! There are only 4 methods available."
    end if

  end subroutine mainchannelRouting

!-----------------------------------------------------------------------

  subroutine Routing_KW(iunit, nt, theDeltaT)
  ! !DESCRIPTION: classic kinematic wave routing method
    implicit none    
    
    integer, intent(in) :: iunit, nt   
    real(r8), intent(in) :: theDeltaT    
    integer  :: k
    real(r8) :: temp_gwl, temp_dwr, temp_gwl0
    character(len=*),parameter :: subname = '(Routing_KW)'

    ! estimate the inflow from upstream units
    TRunoff%erin(iunit,nt) = 0._r8

! tcraig, moved this out of the inner main channel loop to before main channel call
! now it's precomputed as TRunoff%eroutUp
!    do k=1,TUnit%nUp(iunit)
!       TRunoff%erin(iunit,nt) = TRunoff%erin(iunit,nt) - TRunoff%erout(TUnit%iUp(iunit,k),nt)
!    end do
    TRunoff%erin(iunit,nt) = TRunoff%erin(iunit,nt) - TRunoff%eroutUp(iunit,nt)

    ! estimate the outflow
    if(TUnit%rlen(iunit) <= 0._r8) then ! no river network, no channel routing
       TRunoff%vr(iunit,nt) = 0._r8
       TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
    else
       if(TUnit%areaTotal2(iunit)/TUnit%rwidth(iunit)/TUnit%rlen(iunit) > 1e6_r8) then
          TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
       else
!        !TRunoff%vr(iunit,nt) = CRVRMAN(TUnit%rslp(iunit), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
          TRunoff%vr(iunit,nt) = CRVRMAN_nosqrt(TUnit%rslpsqrt(iunit), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
          TRunoff%erout(iunit,nt) = -TRunoff%vr(iunit,nt) * TRunoff%mr(iunit,nt)
          if(-TRunoff%erout(iunit,nt) > TINYVALUE .and. TRunoff%wr(iunit,nt) + (TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT < TINYVALUE) then
             TRunoff%erout(iunit,nt) = -(TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%wr(iunit,nt) / theDeltaT)
             if(TRunoff%mr(iunit,nt) > 0._r8) then
                TRunoff%vr(iunit,nt) = -TRunoff%erout(iunit,nt) / TRunoff%mr(iunit,nt)
             end if
          end if
       end if
    end if

    temp_dwr = TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)
    temp_gwl = TRunoff%qgwl(iunit,nt) * TUnit%area(iunit) * TUnit%frac(iunit)
    temp_gwl0 = temp_gwl
    if(abs(temp_gwl) <= TINYVALUE) then
!       write(iulog,*) 'mosart: ERROR dropping temp_gwl too small'
!       call shr_sys_abort('mosart: ERROR temp_gwl too small')
       temp_gwl = 0._r8
    end if 
    if(temp_gwl < -TINYVALUE) then 
       write(iulog,*) 'mosart: ERROR temp_gwl negative',iunit,nt,TRunoff%qgwl(iunit,nt)
       call shr_sys_abort('mosart: ERROR temp_gwl negative ')
       if(TRunoff%wr(iunit,nt) < TINYVALUE) then
          temp_gwl = 0._r8
       else 
          if(TRunoff%wr(iunit,nt)/theDeltaT + temp_dwr + temp_gwl < -TINYVALUE) then
          !write(iulog,*) 'adjust! ', temp_gwl, -(temp_dwr+TRunoff%wr(iunit,nt)/theDeltaT)
             temp_gwl = -(temp_dwr + TRunoff%wr(iunit,nt) / theDeltaT)
          end if
       end if
    end if
           
    TRunoff%dwr(iunit,nt) = TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt) + temp_gwl

! check for stability
!    if(TRunoff%vr(iunit,nt) < -TINYVALUE .or. TRunoff%vr(iunit,nt) > 30) then
!       write(iulog,*) "Numerical error inRouting_KW, ", iunit,nt,TRunoff%vr(iunit,nt)
!    end if

! check for negative wr
!    if(TRunoff%wr(iunit,nt) > 1._r8 .and. (TRunoff%wr(iunit,nt)/theDeltaT + TRunoff%dwr(iunit,nt))/TRunoff%wr(iunit,nt) < -TINYVALUE) then
!       write(iulog,*) 'negative wr!', TRunoff%wr(iunit,nt), TRunoff%dwr(iunit,nt), temp_dwr, temp_gwl, temp_gwl0, theDeltaT
!       stop          
!    end if 

  end subroutine Routing_KW

!-----------------------------------------------------------------------

  subroutine Routing_MC(iunit, nt, theDeltaT)
  ! !DESCRIPTION: Muskingum-Cunge routing method
    implicit none    
    integer, intent(in) :: iunit, nt   
    real(r8), intent(in) :: theDeltaT
    character(len=*),parameter :: subname = '(Routing_MC)'
   
  end subroutine Routing_MC

!-----------------------------------------------------------------------

  subroutine Routing_THREW(iunit, nt, theDeltaT)
  ! !DESCRIPTION: kinematic wave routing method from THREW model
    implicit none    
    integer, intent(in) :: iunit, nt
    real(r8), intent(in) :: theDeltaT
    character(len=*),parameter :: subname = '(Routing_THREW)'
   
  end subroutine Routing_THREW

!-----------------------------------------------------------------------

  subroutine Routing_DW(iunit, nt, theDeltaT)
  ! !DESCRIPTION: classic diffusion wave routing method
    implicit none    
    integer, intent(in) :: iunit, nt
    real(r8), intent(in) :: theDeltaT
    character(len=*),parameter :: subname = '(Routing_DW)'
   
  end subroutine Routing_DW

!-----------------------------------------------------------------------

  subroutine updateState_hillslope(iunit,nt)
  ! !DESCRIPTION: update the state variables at hillslope
    implicit none    
    integer, intent(in) :: iunit, nt
    character(len=*),parameter :: subname = '(updateState_hillslope)'

    TRunoff%yh(iunit,nt) = TRunoff%wh(iunit,nt) !/ TUnit%area(iunit) / TUnit%frac(iunit) 

  end subroutine updateState_hillslope

!-----------------------------------------------------------------------

  subroutine updateState_subnetwork(iunit,nt)
  ! !DESCRIPTION: update the state variables in subnetwork channel
    implicit none    
    integer, intent(in) :: iunit,nt
    character(len=*),parameter :: subname = '(updateState_subnetwork)'

    if (TUnit%tlen(iunit) > 0._r8 .and. TRunoff%wt(iunit,nt) > 0._r8) then
       TRunoff%mt(iunit,nt) = GRMR(TRunoff%wt(iunit,nt), TUnit%tlen(iunit)) 
       TRunoff%yt(iunit,nt) = GRHT(TRunoff%mt(iunit,nt), TUnit%twidth(iunit))
       TRunoff%pt(iunit,nt) = GRPT(TRunoff%yt(iunit,nt), TUnit%twidth(iunit))
       TRunoff%rt(iunit,nt) = GRRR(TRunoff%mt(iunit,nt), TRunoff%pt(iunit,nt))
    else
       TRunoff%mt(iunit,nt) = 0._r8
       TRunoff%yt(iunit,nt) = 0._r8
       TRunoff%pt(iunit,nt) = 0._r8
       TRunoff%rt(iunit,nt) = 0._r8
    end if

  end subroutine updateState_subnetwork

!-----------------------------------------------------------------------

  subroutine updateState_mainchannel(iunit, nt)
  ! !DESCRIPTION: update the state variables in main channel
    implicit none    
    integer, intent(in) :: iunit, nt
    character(len=*),parameter :: subname = '(updateState_mainchannel)'

    if (TUnit%rlen(iunit) > 0._r8 .and. TRunoff%wr(iunit,nt) > 0._r8) then
       TRunoff%mr(iunit,nt) = GRMR(TRunoff%wr(iunit,nt), TUnit%rlen(iunit)) 
       TRunoff%yr(iunit,nt) = GRHR(TRunoff%mr(iunit,nt), TUnit%rwidth(iunit), TUnit%rwidth0(iunit), TUnit%rdepth(iunit))
       TRunoff%pr(iunit,nt) = GRPR(TRunoff%yr(iunit,nt), TUnit%rwidth(iunit), TUnit%rwidth0(iunit), TUnit%rdepth(iunit))
       TRunoff%rr(iunit,nt) = GRRR(TRunoff%mr(iunit,nt), TRunoff%pr(iunit,nt))
    else
       TRunoff%mr(iunit,nt) = 0._r8
       TRunoff%yr(iunit,nt) = 0._r8
       TRunoff%pr(iunit,nt) = 0._r8
       TRunoff%rr(iunit,nt) = 0._r8
    end if

  end subroutine updateState_mainchannel

!-----------------------------------------------------------------------
    
  function CRVRMAN(slp_, n_, rr_) result(v_)
  ! ! Function for calculating channel velocity according to Manning's equation.
    implicit none
    real(r8), intent(in) :: slp_, n_, rr_ ! slope, manning's roughness coeff., hydraulic radius
    real(r8)             :: v_            ! v_ is  discharge
    
    real(r8) :: ftemp,vtemp
    character(len=*),parameter :: subname = '(CRVRMAN)'

    if(rr_ <= 0._r8) then
       v_ = 0._r8
    else
!tcraig, original code
!       ftemp = 2._r8/3._r8
!       v_ = (rr_**ftemp) * sqrt(slp_) / n_  
!tcraig, produces same answer as original in same time
!       v_ = (rr_**(2._r8/3._r8)) * sqrt(slp_) / n_  

!tcraig, this is faster but NOT bit-for-bit
       v_ = ((rr_*rr_)**(1._r8/3._r8)) * sqrt(slp_) / n_

!debug       if (abs(vtemp - v_)/vtemp > 1.0e-14) then
!debug          write(iulog,*) 'tcx check crvrman ',vtemp, v_
!debug       endif
    end if

  end function CRVRMAN

!-----------------------------------------------------------------------
    
  function CRVRMAN_nosqrt(sqrtslp_, n_, rr_) result(v_)
  ! ! Function for calculating channel velocity according to Manning's equation.
    implicit none
    real(r8), intent(in) :: sqrtslp_, n_, rr_ ! sqrt(slope), manning's roughness coeff., hydraulic radius
    real(r8)             :: v_            ! v_ is  discharge
    
    real(r8) :: ftemp, vtemp
    character(len=*),parameter :: subname = '(CRVRMAN_nosqrt)'

    if(rr_ <= 0._r8) then
       v_ = 0._r8
    else
!tcraig, original code
!       ftemp = 2._r8/3._r8
!       v_ = (rr_**ftemp) * sqrtslp_ / n_  
!tcraig, produces same answer as original in same time
!       v_ = (rr_**(2._r8/3._r8)) * sqrtslp_ / n_  

!tcraig, this is faster but NOT bit-for-bit
       v_ = ((rr_*rr_)**(1._r8/3._r8)) * sqrtslp_ / n_

!debug       if (abs(vtemp - v_)/vtemp > 1.0e-14) then
!debug          write(iulog,*) 'tcx check crvrman_nosqrt ',vtemp, v_
!debug       endif
    end if

  end function CRVRMAN_nosqrt

!-----------------------------------------------------------------------

  function CREHT(hslp_, nh_, Gxr_, yh_) result(eht_)
  ! ! Function for overland from hillslope into the sub-network channels
    implicit none
    real(r8), intent(in) :: hslp_, nh_, Gxr_, yh_ ! topographic slope, manning's roughness coeff., drainage density, overland flow depth
    real(r8)                   :: eht_            ! velocity, specific discharge
    
    real(r8) :: vh_
    character(len=*),parameter :: subname = '(CREHT)'

    vh_ = CRVRMAN(hslp_,nh_,yh_)
    eht_ = Gxr_*yh_*vh_

  end function CREHT

!-----------------------------------------------------------------------

  function CREHT_nosqrt(sqrthslp_, nh_, Gxr_, yh_) result(eht_)
  ! ! Function for overland from hillslope into the sub-network channels
    implicit none
    real(r8), intent(in) :: sqrthslp_, nh_, Gxr_, yh_ ! topographic slope, manning's roughness coeff., drainage density, overland flow depth
    real(r8)                   :: eht_            ! velocity, specific discharge
    
    real(r8) :: vh_
    character(len=*),parameter :: subname = '(CREHT_nosqrt)'

    vh_ = CRVRMAN_nosqrt(sqrthslp_,nh_,yh_)
    eht_ = Gxr_*yh_*vh_

  end function CREHT_nosqrt

!-----------------------------------------------------------------------

  function GRMR(wr_, rlen_) result(mr_)
  ! ! Function for estimate wetted channel area
    implicit none
    real(r8), intent(in) :: wr_, rlen_      ! storage of water, channel length
    real(r8)             :: mr_             ! wetted channel area
    character(len=*),parameter :: subname = '(GRMR)'
    
    mr_ = wr_ / rlen_

  end function GRMR
  
!-----------------------------------------------------------------------

  function GRHT(mt_, twid_) result(ht_)
  ! ! Function for estimating water depth assuming rectangular channel
    implicit none
    real(r8), intent(in) :: mt_, twid_      ! wetted channel area, channel width
    real(r8)             :: ht_             ! water depth
    character(len=*),parameter :: subname = '(GRHT)'
    
    if(mt_ <= TINYVALUE) then
       ht_ = 0._r8
    else
       ht_ = mt_ / twid_
    end if

  end function GRHT

!-----------------------------------------------------------------------

  function GRPT(ht_, twid_) result(pt_)
  ! ! Function for estimating wetted perimeter assuming rectangular channel
    implicit none
    real(r8), intent(in) :: ht_, twid_      ! water depth, channel width
    real(r8)             :: pt_             ! wetted perimeter
    character(len=*),parameter :: subname = '(GRPT)'
    
    if(ht_ <= TINYVALUE) then
       pt_ = 0._r8
    else
       pt_ = twid_ + 2._r8 * ht_
    end if

  end function GRPT

!-----------------------------------------------------------------------

  function GRRR(mr_, pr_) result(rr_)
  ! ! Function for estimating hydraulic radius
    implicit none
    real(r8), intent(in) :: mr_, pr_        ! wetted area and perimeter
    real(r8)             :: rr_             ! hydraulic radius
    character(len=*),parameter :: subname = '(GRRR)'
    
    if(pr_ <= TINYVALUE) then
       rr_ = 0._r8
    else
       rr_ = mr_ / pr_
    end if

  end function GRRR

!-----------------------------------------------------------------------

  function GRHR(mr_, rwidth_, rwidth0_, rdepth_) result(hr_)
  ! ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
  ! ! here assuming the channel cross-section consists of three parts, from bottom to up,
  ! ! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
  ! ! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
  ! ! part 3 is a rectagular with the width rwid0
    implicit none
    real(r8), intent(in) :: mr_, rwidth_, rwidth0_, rdepth_ ! wetted channel area, channel width, flood plain wid, water depth
    real(r8)             :: hr_                             ! water depth
    
    real(r8) :: SLOPE1  ! slope of flood plain, TO DO
    real(r8) :: deltamr_
    character(len=*),parameter :: subname = '(GRHR)'

    SLOPE1 = SLOPE1def
    if(mr_ <= TINYVALUE) then
       hr_ = 0._r8
    else
       if(mr_ - rdepth_*rwidth_ <= TINYVALUE) then ! not flooded
          hr_ = mr_/rwidth_
       else ! if flooded, the find out the equivalent depth
          if(mr_ > rdepth_*rwidth_ + (rwidth_ + rwidth0_)*SLOPE1*((rwidth0_-rwidth_)/2._r8)/2._r8 + TINYVALUE) then
             deltamr_ = mr_ - rdepth_*rwidth_ - (rwidth_ + rwidth0_)*SLOPE1*((rwidth0_ - rwidth_)/2._r8)/2._r8;
             hr_ = rdepth_ + SLOPE1*((rwidth0_ - rwidth_)/2._r8) + deltamr_/(rwidth0_);
          else
             deltamr_ = mr_ - rdepth_*rwidth_;
!           !hr_ = rdepth_ + (-rwidth_+sqrt( rwidth_**2._r8  +4._r8*deltamr_/SLOPE1))*SLOPE1/2._r8
             hr_ = rdepth_ + (-rwidth_+sqrt((rwidth_*rwidth_)+4._r8*deltamr_/SLOPE1))*SLOPE1/2._r8
          end if
       end if
    end if

  end function GRHR
  
!-----------------------------------------------------------------------

  function GRPR(hr_, rwidth_, rwidth0_,rdepth_) result(pr_)
  ! ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
  ! ! here assuming the channel cross-section consists of three parts, from bottom to up,
  ! ! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
  ! ! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
  ! ! part 3 is a rectagular with the width rwid0
    implicit none
    real(r8), intent(in) :: hr_, rwidth_, rwidth0_, rdepth_ ! wwater depth, channel width, flood plain wid, water depth
    real(r8)             :: pr_                             ! water depth
    
    real(r8) :: SLOPE1  ! slope of flood plain, TO DO
    real(r8) :: deltahr_
    logical, save :: first_call = .true.
    character(len=*),parameter :: subname = '(GRPR)'

    SLOPE1 = SLOPE1def
    if (first_call) then
       sinatanSLOPE1defr = 1.0_r8/(sin(atan(SLOPE1def)))
    endif
    first_call = .false.

    if(hr_ < TINYVALUE) then
       pr_ = 0._r8
    else
       if(hr_ <= rdepth_ + TINYVALUE) then ! not flooded
          pr_ = rwidth_ + 2._r8*hr_
       else
          if(hr_ > rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1 + TINYVALUE) then
             deltahr_ = hr_ - rdepth_ - ((rwidth0_-rwidth_)/2._r8)*SLOPE1
!           !pr_ = rwidth_ + 2._r8*(rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1/sin(atan(SLOPE1)) + deltahr_)
             pr_ = rwidth_ + 2._r8*(rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1*sinatanSLOPE1defr + deltahr_)
          else
!           !pr_ = rwidth_ + 2._r8*(rdepth_ + (hr_ - rdepth_)/sin(atan(SLOPE1)))
             pr_ = rwidth_ + 2._r8*(rdepth_ + (hr_ - rdepth_)*sinatanSLOPE1defr)
          end if
       end if
    end if

  end function GRPR 
  
!-----------------------------------------------------------------------

  subroutine createFile(nio, fname)
  ! !DESCRIPTION: create a new file. if a file with the same name exists, delete it then create a new one
    implicit none
    character(len=*), intent(in) :: fname ! file name
    integer, intent(in) :: nio            !unit of the file to create

    integer :: ios
    logical :: filefound
    character(len=1000) :: cmd
    character(len=*),parameter :: subname = '(createFile)'

    inquire (file=fname, exist=filefound)
    if(filefound) then
       cmd = 'rm '//trim(fname)
       call system(cmd)
    end if
    open (unit=nio, file=fname, status="new", action="write", iostat=ios)
    if(ios /= 0) then
       print*, "cannot create file ", fname
    end if

  end subroutine createFile
  
!-----------------------------------------------------------------------

  subroutine printTest(nio)
  ! !DESCRIPTION: output the simulation results into external files
    implicit none
    integer, intent(in) :: nio        ! unit of the file to print
    
    integer :: IDlist(1:5) = (/151,537,687,315,2080/)
    integer :: nt
    integer :: ios,ii                    ! flag of io status
    character(len=*),parameter :: subname = '(printTest)'

    write(unit=nio,fmt="(15(e20.11))") TRunoff%etin(IDlist(1),1)/TUnit%area(IDlist(1)), TRunoff%erlateral(IDlist(1),1)/TUnit%area(IDlist(1)), TRunoff%flow(IDlist(1),1), &
                                       TRunoff%etin(IDlist(2),1)/TUnit%area(IDlist(2)), TRunoff%erlateral(IDlist(2),1)/TUnit%area(IDlist(2)), TRunoff%flow(IDlist(2),1), &
                     TRunoff%etin(IDlist(3),1)/TUnit%area(IDlist(3)), TRunoff%erlateral(IDlist(3),1)/TUnit%area(IDlist(3)), TRunoff%flow(IDlist(3),1), &
                     TRunoff%etin(IDlist(4),1)/TUnit%area(IDlist(4)), TRunoff%erlateral(IDlist(4),1)/TUnit%area(IDlist(4)), TRunoff%flow(IDlist(4),1), &
                     TRunoff%etin(IDlist(5),1)/TUnit%area(IDlist(5)), TRunoff%erlateral(IDlist(5),1)/TUnit%area(IDlist(5)), TRunoff%flow(IDlist(5),1)
    !nt = 1
    !write(unit=nio,fmt="((a10),(e20.11))") theTime, TRunoff%flow(ii,nt)
    !write(unit=nio,fmt="((a10),6(e20.11))") theTime, TRunoff%qsur(ii,nt), TRunoff%qsub(ii,nt), TRunoff%etin(ii,nt)/(TUnit%area(ii)*TUnit%frac(ii)), TRunoff%erlateral(ii,nt)/(TUnit%area(ii)*TUnit%frac(ii)), TRunoff%erin(ii,nt), TRunoff%flow(ii,nt)
    !if(TRunoff%yr(ii,nt) > 0._r8) then
    !    write(unit=nio,fmt="((a10),6(e20.11))") theTime, TRunoff%mr(ii,nt)/TRunoff%yr(ii,nt),TRunoff%yr(ii,nt), TRunoff%vr(ii,nt), TRunoff%erin(ii,nt), TRunoff%erout(ii,nt)/(TUnit%area(ii)*TUnit%frac(ii)), TRunoff%flow(ii,nt)
      !else
    !    write(unit=nio,fmt="((a10),6(e20.11))") theTime, TRunoff%mr(ii,nt)-TRunoff%mr(ii,nt),TRunoff%yr(ii,nt), TRunoff%vr(ii,nt), TRunoff%erin(ii,nt), TRunoff%erout(ii,nt)/(TUnit%area(ii)*TUnit%frac(ii)), TRunoff%flow(ii,nt)
    !end if
    !write(unit=nio,fmt="((a10),7(e20.11))") theTime, TRunoff%erlateral(ii,nt)/(TUnit%area(ii)*TUnit%frac(ii)), TRunoff%wr(ii,nt),TRunoff%mr(ii,nt), TRunoff%yr(ii,nt), TRunoff%pr(ii,nt), TRunoff%rr(ii,nt), TRunoff%flow(ii,nt)
    !write(unit=nio,fmt="((a10),7(e20.11))") theTime, TRunoff%yh(ii,nt), TRunoff%dwh(ii,nt),TRunoff%etin(ii,nt), TRunoff%vr(ii,nt), TRunoff%erin(ii,nt), TRunoff%erout(ii,nt)/(TUnit%area(ii)*TUnit%frac(ii)), TRunoff%flow(ii,nt)
  
  end subroutine printTest

!-----------------------------------------------------------------------

end MODULE MOSART_physics_mod

