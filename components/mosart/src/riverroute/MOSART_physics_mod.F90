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
  use RtmVar        , only : iulog, barrier_timers, wrmflag, inundflag, sediflag, heatflag, rstraflag, use_ocn_rof_two_way
  use RunoffMod     , only : Tctl, TUnit, TRunoff, Theat, TPara, rtmCTL, &
                             SMatP_upstrm, avsrc_upstrm, avdst_upstrm, SMatP_dnstrm, avsrc_dnstrm, avdst_dnstrm
  use MOSART_heat_mod
  use MOSART_stra_mod
  use RtmSpmd       , only : masterproc, mpicom_rof, iam
  use RtmTimeManager, only : get_curr_date, is_new_month

  use WRM_type_mod  , only : ctlSubwWRM, WRMUnit, StorWater
  use WRM_modules   , only : irrigationExtractionSubNetwork, &
                             irrigationExtractionMainChannel, &
                             Regulation, ExtractionRegulatedFlow
  use WRM_returnflow, only : insert_returnflow_channel, &
                             insert_returnflow_soilcolumn, &
                             estimate_returnflow_deficit
  use WRM_subw_io_mod, only : WRM_readDemand, WRM_computeRelease
  use MOSARTinund_Core_MOD, only: ChnlFPexchg
  use rof_cpl_indices, only : nt_rtm, rtm_tracers, nt_nliq, nt_nice, nt_nmud, nt_nsan, KW, DW
  use perf_mod, only: t_startf, t_stopf
  use mct_mod
  use MOSART_BGC_type, only : TSedi
  use MOSART_sediment_mod
  use MOSART_RES_type
  use MOSART_reservoir_mod

  implicit none
  private

  real(r8), parameter :: TINYVALUE = 1.0e-14_r8  ! double precision variable has a significance of about 16 decimal digits
  integer  :: nt               ! loop indices
  real(r8), parameter :: SLOPE1def = 0.1_r8        ! here give it a small value in order to avoid the abrupt change of hydraulic radidus etc.
  real(r8) :: sinatanSLOPE1defr   ! 1.0/sin(atan(slope1))
  real(r8), parameter :: MaxStorageDepleted = 0.95_r8        ! maximum storage allowed to deplete in a single step -- a trick to keep water balance and numerical stability
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
    
    integer :: iunit, idam, m, k, unitUp, cnt, ier, dd, nSubStep   !local index
    real(r8) :: temp_erout, localDeltaT, temp_haout, temp_Tt, temp_Tr, temp_T, temp_ha
    real(r8) :: mud_erout, san_erout, temp_ehexch, temp_etexch, temp_erexch
    real(r8) :: negchan 
    integer  :: numSubSteps
    integer  :: yr,mon,day,tod
    real(r8) :: myTINYVALUE
    character(len=*),parameter :: subname = '(Euler)'
    real(r8) :: tmp1, tmp2
    !------------------

    myTINYVALUE = 1.e-6

    call get_curr_date(yr, mon, day, tod)
    !------------------
    ! WRM prep
    !------------------


    if (wrmflag) then
       if ( ctlSubwWRM%ReturnFlowFlag > 0) then
          call insert_returnflow_soilcolumn
       endif
       !call readPotentialEvap(trim(theTime))
       if ( is_new_month() ) then
         if (masterproc) write(iulog,*) trim(subname),' updating monthly data at ',yr,mon,day,tod
         if (ctlSubwWRM%ExternalDemandFlag > 0) then
          call WRM_readDemand()
         end if 
         call WRM_computeRelease() !! about regulation
       end if

        if (ctlSubwWRM%ExtractionFlag > 0) then 
         do iunit=rtmCTL%begr,rtmCTL%endr
           if (ctlSubwWRM%ExternalDemandFlag == 0) then  ! if demand is from ELM, reset the demand0 every timestep
             StorWater%demand0(iunit) = 0
           endif
            do nt=nt_nliq,nt_nice
              if (TUnit%mask(iunit) > 0) then
                  if (ctlSubwWRM%ExternalDemandFlag == 0) then  ! if demand is from ELM
                      StorWater%demand0(iunit) = StorWater%demand0(iunit) - TRunoff%qdem(iunit,nt) * TUnit%area(iunit) * TUnit%frac(iunit)
                  endif
              endif
            enddo  
          enddo
       end if
       StorWater%demand = StorWater%demand0 * Tctl%DeltaT
       !supply is set to zero in RtmMod so it can be accumulated there for the budget
       !StorWater%supply = 0._r8
       !StorWater%deficit =0._r8
    endif

    !------------------
    ! hillslope
    !------------------

    call t_startf('mosartr_hillslope')
    do nt=nt_nliq,nt_nice
    if (TUnit%euler_calc(nt)) then
    do iunit=rtmCTL%begr,rtmCTL%endr
       if(TUnit%mask(iunit) > 0) then
          call hillslopeRouting(iunit,nt,Tctl%DeltaT)
          TRunoff%wh(iunit,nt) = TRunoff%wh(iunit,nt) + TRunoff%dwh(iunit,nt) * Tctl%DeltaT
          call UpdateState_hillslope(iunit,nt)
          TRunoff%etin(iunit,nt) = (-TRunoff%ehout(iunit,nt) + TRunoff%qsub(iunit,nt)) * TUnit%area(iunit) * TUnit%frac(iunit)
          if (heatflag) then
          if (nt==nt_nliq) then
              call hillslopeHeat(iunit, Tctl%DeltaT)
          end if
          end if
       endif
    end do
    endif
    end do

    if (sediflag .and. TUnit%euler_calc(nt_nmud)) then
    do iunit=rtmCTL%begr,rtmCTL%endr
       if(TUnit%mask(iunit) > 0) then
          call hillslopeSediment(iunit, Tctl%DeltaT)
          TRunoff%etin(iunit,nt_nmud) = (-TRunoff%ehout(iunit,nt_nmud) + TRunoff%qsub(iunit,nt_nmud)) * TUnit%area(iunit) * TUnit%frac(iunit)
          TRunoff%ehexch_avg(iunit,nt_nmud) = 0._r8

          !! note: only when the soil erosion subroutine is turned on, the ehexchange item is meaningful, otherwise always zero
          !call soilErosion(iunit, Tctl%DeltaT)
          !TRunoff%ehexch_avg(iunit,nt_nmud) = TRunoff%etin(iunit,nt_nmud)
       endif
    end do
    endif
    call t_stopf('mosartr_hillslope')

    TRunoff%flow = 0._r8
    TRunoff%erowm_regi = 0._r8
    TRunoff%erowm_regf = 0._r8
    TRunoff%eroup_lagi = 0._r8
    TRunoff%eroup_lagf = 0._r8
    TRunoff%eroutup_avg = 0._r8
    TRunoff%erlat_avg = 0._r8
    if (heatflag) then
       THeat%Ha_eroutup_avg = 0._r8
       THeat%Ha_erlat_avg = 0._r8
       THeat%Tt_avg = 0._r8
       THeat%Tr_avg = 0._r8
    endif
    if (inundflag) then
       TRunoff%se_rf = 0._r8
    endif

    TRunoff%etexch_avg = 0._r8
    TRunoff%erexch_avg = 0._r8
    negchan = 9999.0_r8

    ! subcycling within MOSART begins
    do m=1,Tctl%DLevelH2R

       !------------------
       ! subnetwork
       !------------------

       call t_startf('mosartr_subnetwork')    
       TRunoff%erlateral(:,:) = 0._r8
       if (heatflag) THeat%ha_lateral(:) = 0._r8
       TRunoff%etexchange = 0._r8
       do nt=nt_nliq,nt_nice ! water transport
       if (TUnit%euler_calc(nt)) then
       do iunit=rtmCTL%begr,rtmCTL%endr
          temp_Tt = 0._r8
          if(TUnit%mask(iunit) > 0) then
!extraction from subnetwork here from wt
             if (wrmflag) then
                if (nt == nt_nliq) then
                   if  (ctlSubwWRM%ExtractionFlag > 0 .and. TRunoff%yt(iunit,nt_nliq) >= 0.1_r8) then
                      localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
                      call irrigationExtractionSubNetwork(iunit, localDeltaT )
                      call UpdateState_subnetwork(iunit,nt)
                   endif
                endif
             endif
             localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_t(iunit)
             do k=1,TUnit%numDT_t(iunit)
                call subnetworkRouting(iunit,nt,localDeltaT)
                TRunoff%wt(iunit,nt) = TRunoff%wt(iunit,nt) + TRunoff%dwt(iunit,nt) * localDeltaT
                call UpdateState_subnetwork(iunit,nt)
                TRunoff%erlateral(iunit,nt) = TRunoff%erlateral(iunit,nt)-TRunoff%etout(iunit,nt)
                if (heatflag) then
                if (nt==nt_nliq) then
                    if(TUnit%tlen(iunit) > myTINYVALUE) then
                          if(TRunoff%yt(iunit,nt_nliq) >= 0.2_r8) then 
                              call subnetworkHeat(iunit,localDeltaT)
                              call subnetworkTemp(iunit)
                          elseif(TRunoff%yt(iunit,nt_nliq) <= 0.05_r8) then
                              call subnetworkHeat_simple(iunit,localDeltaT)
                              THeat%Tt(iunit) = cr_S_curve(iunit,THeat%forc_t(iunit))
                          else
                              temp_T = 0._r8
                              temp_ha = 0._r8
                              nSubStep = 10
                              do dd=1,nSubStep
                                  call subnetworkHeat(iunit,localDeltaT/nSubStep)
                                  call subnetworkTemp(iunit)
                                  temp_T = temp_T + THeat%Tt(iunit)
                                  temp_ha = temp_ha + THeat%Ha_t2r(iunit)
                              end do
                              THeat%Tt(iunit) = temp_T/nSubStep
                              THeat%Ha_t2r(iunit) = temp_ha/nSubStep
                          end if
                          THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) - THeat%Ha_t2r(iunit)
                          temp_Tt = temp_Tt + THeat%Tt(iunit)
                      else
                          call subnetworkHeat_simple(iunit,localDeltaT)
                          call subnetworkTemp_simple(iunit)
                          THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) - THeat%Ha_t2r(iunit)
                          temp_Tt = temp_Tt + THeat%Tt(iunit)
                    end if
                end if
                end if
             end do ! numDT_t
             TRunoff%erlateral(iunit,nt) = TRunoff%erlateral(iunit,nt) / TUnit%numDT_t(iunit)
             if (heatflag) then
             if (nt==nt_nliq) then
                 THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) / TUnit%numDT_t(iunit)
                 temp_Tt = temp_Tt / TUnit%numDT_t(iunit)
                 THeat%Tt_avg(iunit) = THeat%Tt_avg(iunit) + temp_Tt
             end if
             end if
          endif
       end do ! iunit
       endif  ! euler_calc
       end do ! nt

       !! the treatment of mud and san is special since these two are interacting with each other
       !do nt=nmud,nt_nsan ! sediment transport
       if (sediflag .and. TUnit%euler_calc(nt_nmud)) then
       do iunit=rtmCTL%begr,rtmCTL%endr
          if(TUnit%mask(iunit) > 0) then
             localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_t(iunit)
             do k=1,TUnit%numDT_t(iunit)
                call subnetworkSediment(iunit,localDeltaT)
                TRunoff%wt(iunit,nt_nmud) = TRunoff%wt(iunit,nt_nmud) + TRunoff%dwt(iunit,nt_nmud) * localDeltaT
                TRunoff%wt(iunit,nt_nsan) = TRunoff%wt(iunit,nt_nsan) + TRunoff%dwt(iunit,nt_nsan) * localDeltaT
                TRunoff%wt_al(iunit,nt_nmud) = TRunoff%wt_al(iunit,nt_nmud) + TRunoff%dwt_al(iunit,nt_nmud) * localDeltaT
                TRunoff%wt_al(iunit,nt_nsan) = TRunoff%wt_al(iunit,nt_nsan) + TRunoff%dwt_al(iunit,nt_nsan) * localDeltaT
                call UpdateState_subnetwork(iunit,nt_nmud)
                call UpdateState_subnetwork(iunit,nt_nsan)
                TRunoff%erlateral(iunit,nt_nmud) = TRunoff%erlateral(iunit,nt_nmud) - TRunoff%etout(iunit,nt_nmud)
                TRunoff%erlateral(iunit,nt_nsan) = TRunoff%erlateral(iunit,nt_nsan) - TRunoff%etout(iunit,nt_nsan)
                TRunoff%etexchange(iunit,nt_nmud) = TRunoff%etexchange(iunit,nt_nmud) + TSedi%ermb_t(iunit)
                TRunoff%etexchange(iunit,nt_nsan) = TRunoff%etexchange(iunit,nt_nsan) + TSedi%ersb_t(iunit)
             end do ! numDT_t
             TRunoff%erlateral(iunit,nt_nmud) = TRunoff%erlateral(iunit,nt_nmud) / TUnit%numDT_t(iunit)
             TRunoff%erlateral(iunit,nt_nsan) = TRunoff%erlateral(iunit,nt_nsan) / TUnit%numDT_t(iunit)
             TRunoff%etexchange(iunit,nt_nmud) = TRunoff%etexchange(iunit,nt_nmud) / TUnit%numDT_t(iunit)
             TRunoff%etexchange(iunit,nt_nsan) = TRunoff%etexchange(iunit,nt_nsan) / TUnit%numDT_t(iunit)
             TRunoff%etexch_avg(iunit,nt_nmud) = TRunoff%etexch_avg(iunit,nt_nmud) + TRunoff%etexchange(iunit,nt_nmud)
             TRunoff%etexch_avg(iunit,nt_nsan) = TRunoff%etexch_avg(iunit,nt_nsan) + TRunoff%etexchange(iunit,nt_nsan)
          endif

!#ifdef INCLUDE_WRM
          !! TODO: sediment trapping by small reservoirs on sub-network channel
          !if (sediflag .and. wrmflag) then
          !   localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
          !   do nt=nt_nmud,nt_nsan  ! I need to use something similar to storwater%supply and storwater%demand to keep the budget right
          !      !TRunoff%erowm_regi(iunit,nt) = TRunoff%erowm_regi(iunit,nt) + TRunoff%erlateral(iunit,nt)
          !   enddo
          !      
          !   call res_trapping_t(iunit,nt_nmud)
          !   Tres%wres_t(iunit,nt_nmud) = Tres%wres_t(iunit,nt_nmud) + Tres%dwres_t(iunit,nt_nmud) * localDeltaT
          !   call res_trapping_t(iunit,nt_nsan)
          !   Tres%wres_t(iunit,nt_nsan) = Tres%wres_t(iunit,nt_nsan) + Tres%dwres_t(iunit,nt_nsan) * localDeltaT
          !    
          !   do nt=nt_nmud,nt_nsan
          !     !TRunoff%erowm_regf(iunit,nt) = TRunoff%erowm_regf(iunit,nt) + TRunoff%erlateral(iunit,nt)
          !   enddo
          !end if                 
!#endif
       end do ! iunit
       endif  ! euler_calc

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
       if (heatflag) THeat%Ha_eroutUp = 0._r8
#ifdef NO_MCT
       do iunit=rtmCTL%begr,rtmCTL%endr
       do k=1,rtmCTL%nUp(iunit)
          unitUp = rtmCTL%iUp(iunit,k)
          do nt=1,nt_rtm
             TRunoff%eroutUp(iunit,nt) = TRunoff%eroutUp(iunit,nt) + TRunoff%erout(unitUp,nt)
          end do
          if (heatflag) then
              THeat%Ha_eroutUp(iunit) = THeat%Ha_eroutUp(iunit) + THeat%Ha_rout(unitUp)
          end if
       end do
       end do
#else
       !--- copy erout into avsrc_upstrm ---
       call mct_avect_zero(avsrc_upstrm)
       cnt = 0
       do iunit = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          do nt = 1,nt_rtm
             avsrc_upstrm%rAttr(nt,cnt) = TRunoff%erout(iunit,nt)
          enddo
          if (heatflag) then
              avsrc_upstrm%rAttr(nt_rtm+1,cnt) = THeat%Ha_rout(iunit)
          end if
       enddo
       call mct_avect_zero(avdst_upstrm)

       call mct_sMat_avMult(avsrc_upstrm, sMatP_upstrm, avdst_upstrm)

       !--- add mapped eroutUp to TRunoff ---
       cnt = 0
       do iunit = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          do nt = 1,nt_rtm
             TRunoff%eroutUp(iunit,nt) = avdst_upstrm%rAttr(nt,cnt)
          enddo
          if (heatflag) then
              THeat%Ha_eroutUp(iunit) = avdst_upstrm%rAttr(nt_rtm+1,cnt)
          end if
       enddo
       
       if (Tctl%RoutingMethod == DW ) then
          ! retrieve water depth in downstream channels
          call mct_aVect_zero(avsrc_dnstrm)
          cnt = 0
          do iunit = rtmCTL%begr,rtmCTL%endr
             cnt = cnt + 1
             avsrc_dnstrm%rAttr(nt_nliq,cnt) = TRunoff%yr(iunit,nt_nliq)
          enddo
          call mct_aVect_zero(avdst_dnstrm)
          call mct_sMat_avMult(avsrc_dnstrm, sMatP_dnstrm, avdst_dnstrm)
          cnt = 0
          do iunit = rtmCTL%begr,rtmCTL%endr
             cnt = cnt + 1
             TRunoff%yr_dstrm(iunit) = avdst_dnstrm%rAttr(nt_nliq,cnt)
          enddo
          
          ! retrieve water storage in downstream channels
          call mct_aVect_zero(avsrc_dnstrm)
          cnt = 0
          do iunit = rtmCTL%begr,rtmCTL%endr
             cnt = cnt + 1
             do nt = 1,nt_rtm
                 avsrc_dnstrm%rAttr(nt,cnt) = TRunoff%wr(iunit,nt)
             enddo
          enddo
          call mct_aVect_zero(avdst_dnstrm)
          call mct_sMat_avMult(avsrc_dnstrm, sMatP_dnstrm, avdst_dnstrm)
          cnt = 0
          do iunit = rtmCTL%begr,rtmCTL%endr
             cnt = cnt + 1
             do nt = 1,nt_rtm
                 TRunoff%wr_dstrm(iunit,nt) = avdst_dnstrm%rAttr(nt,cnt)
             enddo
          enddo

          ! retrieve concentration of BGC in downstream channels
          call mct_aVect_zero(avsrc_dnstrm)
          cnt = 0
          do iunit = rtmCTL%begr,rtmCTL%endr
             cnt = cnt + 1
             do nt = 1,nt_rtm
                 avsrc_dnstrm%rAttr(nt,cnt) = TRunoff%conc_r(iunit,nt)
             end do
          enddo
          call mct_aVect_zero(avdst_dnstrm)
          call mct_sMat_avMult(avsrc_dnstrm, sMatP_dnstrm, avdst_dnstrm)
          cnt = 0
          do iunit = rtmCTL%begr,rtmCTL%endr
             cnt = cnt + 1
             do nt = 1,nt_rtm
                 TRunoff%conc_r_dstrm(iunit, nt) = avdst_dnstrm%rAttr(nt,cnt)
             end do
          enddo

          ! retrieve total inflow in downstream channels
          call mct_aVect_zero(avsrc_dnstrm)
          cnt = 0
          do iunit = rtmCTL%begr,rtmCTL%endr
             cnt = cnt + 1
             do nt = 1,nt_rtm
                 avsrc_dnstrm%rAttr(nt,cnt) = TRunoff%erin(iunit,nt)
             end do
          enddo
          call mct_aVect_zero(avdst_dnstrm)
          call mct_sMat_avMult(avsrc_dnstrm, sMatP_dnstrm, avdst_dnstrm)
          cnt = 0
          do iunit = rtmCTL%begr,rtmCTL%endr
             cnt = cnt + 1
             do nt = 1,nt_rtm
                 TRunoff%erin_dstrm(iunit, nt) = avdst_dnstrm%rAttr(nt,cnt)
             end do
          enddo

          ! add dstrm BC for rof ocn coupling
          if (use_ocn_rof_two_way) then
             do iunit=rtmCTL%begr,rtmCTL%endr
                if ( (rtmCTL%mask(iunit) .eq. 3) .and. (TUnit%ocn_rof_coupling_ID(iunit) .eq. 1) ) then
                   TRunoff%yr_dstrm(iunit) = rtmCTL%ssh(iunit) + Tunit%vdatum_conversion(iunit) ! assign ocn's water depth to the dstrm component of the specified outlet cell 
                end if
             end do
          end if

       end if
#endif
       call t_stopf('mosartr_SMeroutUp')    

       TRunoff%eroutup_avg = TRunoff%eroutup_avg + TRunoff%eroutUp
       TRunoff%erlat_avg   = TRunoff%erlat_avg   + TRunoff%erlateral
       if (heatflag) then
           THeat%Ha_eroutup_avg = THeat%Ha_eroutup_avg + THeat%Ha_eroutUp
           THeat%Ha_erlat_avg   = THeat%Ha_erlat_avg   + THeat%Ha_lateral
       end if

       !------------------
       ! channel routing
       !------------------

       call t_startf('mosartr_chanroute')    
       TRunoff%erexchange = 0._r8       
       do nt=nt_nliq,nt_nice ! water transport
       if (TUnit%euler_calc(nt)) then
       do iunit=rtmCTL%begr,rtmCTL%endr
          if(TUnit%mask(iunit) > 0) then
             temp_erout = 0._r8
             temp_haout = 0._r8
             temp_Tr = 0._r8
             if(Tctl%RoutingMethod == KW) then  ! local stepping method only applicable for kinamatic wave routing method
                 numSubSteps = TUnit%numDT_r(iunit)
                 localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/numSubSteps
                 do k=1,numSubSteps
                    call mainchannelRouting(iunit,nt,localDeltaT)    
                    TRunoff%wr(iunit,nt) = TRunoff%wr(iunit,nt) + TRunoff%dwr(iunit,nt) * localDeltaT
                    !! check for negative channel storage
                    !if(TRunoff%wr(iunit,1) < -1.e-10) then
                    !   write(iulog,*) 'Negative channel storage! ', iunit, TRunoff%wr(iunit,1), TRunoff%erin(iunit,1), TRunoff%erout(iunit,1), rtmCTL%nUp(iunit)
                    !   call shr_sys_abort('mosart: negative channel storage')
                    !end if
                    call UpdateState_mainchannel(iunit,nt)
                    temp_erout = temp_erout + TRunoff%erout(iunit,nt) ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
                 end do
             elseif(Tctl%RoutingMethod == DW) then ! diffusion wave routing method
                 numSubSteps = 20 ! now set as 20, could be adjusted as needed.
                 localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/numSubSteps
                 TRunoff%rslp_energy(iunit) = CRRSLP(iunit)
                 do k=1,numSubSteps
                    call Leapfrog(iunit,nt,localDeltaT)  ! note updating wr and other states are done in Leapfrog
                    !! check for negative channel storage
                    !if(TRunoff%wr(iunit,1) < -1.e-10) then
                    !   write(iulog,*) 'Negative channel storage! ', iunit, TRunoff%wr(iunit,1), TRunoff%erin(iunit,1), TRunoff%erout(iunit,1), rtmCTL%nUp(iunit)
                    !   call shr_sys_abort('mosart: negative channel storage')
                    !end if
                    temp_erout = temp_erout + TRunoff%erout(iunit,nt) ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
                 end do
             end if
             
             temp_erout = temp_erout / numSubSteps
             TRunoff%erout(iunit,nt) = temp_erout
             if (heatflag) then
             if (nt==nt_nliq) then
                 do k=1,TUnit%numDT_r(iunit)                
                    if(TUnit%rlen(iunit) > myTINYVALUE) then
                        if(TRunoff%yr(iunit,nt_nliq) >= 0.2_r8) then
                            call mainchannelHeat(iunit, localDeltaT)
                            call mainchannelTemp(iunit)
                        elseif(TRunoff%yr(iunit,nt_nliq) <= 0.05_r8) then
                            call mainchannelHeat_simple(iunit, localDeltaT)
                            THeat%Tr(iunit) = cr_S_curve(iunit,THeat%forc_t(iunit))
                        else
                            temp_T = 0._r8
                            temp_ha = 0._r8
                            nSubStep = 10
                            do dd=1,nSubStep
                                call mainchannelHeat(iunit, localDeltaT/nSubStep)
                                call mainchannelTemp(iunit)
                                temp_T = temp_T + THeat%Tr(iunit)
                                temp_ha = temp_ha + THeat%ha_rout(iunit)
                            end do
                            THeat%Tr(iunit) = temp_T/nSubStep
                            THeat%ha_rout(iunit) = temp_ha/nSubStep
                        end if
                        temp_haout = temp_haout + THeat%ha_rout(iunit)
                        temp_Tr = temp_Tr + THeat%Tr(iunit)
                    else
                        call mainchannelHeat_simple(iunit, localDeltaT)
                        call mainchannelTemp_simple(iunit)
                        temp_haout = temp_haout + THeat%ha_rout(iunit)
                        temp_Tr = temp_Tr + THeat%Tr(iunit)
                    end if                
                 end do
                 temp_haout = temp_haout / TUnit%numDT_r(iunit)
                 THeat%ha_rout(iunit) = temp_haout
                 temp_Tr = temp_Tr / TUnit%numDT_r(iunit)
                 THeat%Tr_avg(iunit) = THeat%Tr_avg(iunit) + temp_Tr
             end if
             end if
!#ifdef INCLUDE_WRM
             if (wrmflag) then
                if (nt == nt_nliq) then
                   localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
                   if (ctlSubwWRM%ExtractionMainChannelFlag > 0 .AND. ctlSubwWRM%ExtractionFlag > 0  .and. TRunoff%yr(iunit,nt_nliq) >= 0.1_r8) then
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
                ! moved out of loop
                   if ( ctlSubwWRM%RegulationFlag>0 ) then
                      call Regulation(iunit, localDeltaT)
                      if (heatflag .and. rstraflag) then
                          call stratification(iunit, localDeltaT,nt)
                          call reservoirHeat(iunit, localDeltaT)
                      elseif (heatflag .and. (.not.rstraflag)) then
                          call reservoirHeat(iunit, localDeltaT)
                      end if
                   endif
                endif
                ! do not update wr after regulation or extraction from reservoir release. Because of the regulation, 
                ! the wr might get to crazy uncontrolled values, assume in this case wr is not changed. The storage in reservoir handles it.
             endif

             Trunoff%eroup_lagf(iunit,nt) = Trunoff%eroup_lagf(iunit,nt) - Trunoff%erout(iunit,nt)
             TRunoff%flow(iunit,nt) = TRunoff%flow(iunit,nt) - TRunoff%erout(iunit,nt)
          endif

       end do ! iunit
       endif  ! euler_calc
       end do ! nt      
       !! the mud and sand processes are treated together
       !do nt=nmud,nt_nsan ! sediment transport
       if (sediflag .and. TUnit%euler_calc(nt_nmud)) then
       do iunit=rtmCTL%begr,rtmCTL%endr
          if(TUnit%mask(iunit) > 0) then
             mud_erout = 0._r8
             san_erout = 0._r8
             if(Tctl%RoutingMethod==KW) then  ! local stepping method only applies for the kinematic wave routing method
                 numSubSteps = TUnit%numDT_r(iunit)
                 localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/numSubSteps
                 do k=1,numSubSteps
                    call mainchannelSediment(iunit,localDeltaT)
                    TRunoff%wr(iunit,nt_nmud) = TRunoff%wr(iunit,nt_nmud) + TRunoff%dwr(iunit,nt_nmud) * localDeltaT
                    TRunoff%wr(iunit,nt_nsan) = TRunoff%wr(iunit,nt_nsan) + TRunoff%dwr(iunit,nt_nsan) * localDeltaT
                    TRunoff%wr_al(iunit,nt_nmud) = TRunoff%wr_al(iunit,nt_nmud) + TRunoff%dwr_al(iunit,nt_nmud) * localDeltaT
                    TRunoff%wr_al(iunit,nt_nsan) = TRunoff%wr_al(iunit,nt_nsan) + TRunoff%dwr_al(iunit,nt_nsan) * localDeltaT
                    call UpdateState_mainchannel(iunit,nt_nmud)
                    call UpdateState_mainchannel(iunit,nt_nsan)
                    mud_erout = mud_erout + TRunoff%erout(iunit,nt_nmud) ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
                    san_erout = san_erout + TRunoff%erout(iunit,nt_nsan) ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
                    TRunoff%erexchange(iunit,nt_nmud) = TRunoff%erexchange(iunit,nt_nmud) + TSedi%ermb_r(iunit)
                    TRunoff%erexchange(iunit,nt_nsan) = TRunoff%erexchange(iunit,nt_nsan) + TSedi%ersb_r(iunit)
                 end do
             elseif(Tctl%RoutingMethod==DW) then
                 numSubSteps = 20
                 localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/numSubSteps
                 do k=1,numSubSteps
                    call Leapfrog_sed(iunit,localDeltaT)  ! note updating wr and other states are done in Leapfrog already
                    mud_erout = mud_erout + TRunoff%erout(iunit,nt_nmud) ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
                    san_erout = san_erout + TRunoff%erout(iunit,nt_nsan) ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
                    TRunoff%erexchange(iunit,nt_nmud) = TRunoff%erexchange(iunit,nt_nmud) + TSedi%ermb_r(iunit)
                    TRunoff%erexchange(iunit,nt_nsan) = TRunoff%erexchange(iunit,nt_nsan) + TSedi%ersb_r(iunit)
                 end do
             end if
             mud_erout = mud_erout / numSubSteps
             TRunoff%erout(iunit,nt_nmud) = mud_erout
             san_erout = san_erout / numSubSteps
             TRunoff%erout(iunit,nt_nsan) = san_erout

             Trunoff%eroup_lagf(iunit,nt_nmud) = Trunoff%eroup_lagf(iunit,nt_nmud) - Trunoff%erout(iunit,nt_nmud)
             TRunoff%flow(iunit,nt_nmud) = TRunoff%flow(iunit,nt_nmud) - TRunoff%erout(iunit,nt_nmud)
             TRunoff%erexchange(iunit,nt_nmud) = TRunoff%erexchange(iunit,nt_nmud) / numSubSteps
             TRunoff%erexch_avg(iunit,nt_nmud) = TRunoff%erexch_avg(iunit,nt_nmud) + TRunoff%erexchange(iunit,nt_nmud)

             Trunoff%eroup_lagf(iunit,nt_nsan) = Trunoff%eroup_lagf(iunit,nt_nsan) - Trunoff%erout(iunit,nt_nsan)
             TRunoff%flow(iunit,nt_nsan) = TRunoff%flow(iunit,nt_nsan) - TRunoff%erout(iunit,nt_nsan)
             TRunoff%erexchange(iunit,nt_nsan) = TRunoff%erexchange(iunit,nt_nsan) / numSubSteps
             TRunoff%erexch_avg(iunit,nt_nsan) = TRunoff%erexch_avg(iunit,nt_nsan) + TRunoff%erexchange(iunit,nt_nsan)

          endif

!#ifdef INCLUDE_WRM
          !! Assume that reservoir regulation will only affect suspended load by changing the flow conditions, but do not directly affect sediment flux or storage
          if (sediflag .and. wrmflag) then
             localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
             do nt=nt_nmud,nt_nsan
                TRunoff%erowm_regi(iunit,nt) = TRunoff%erowm_regi(iunit,nt) - TRunoff%erout(iunit,nt)
                TRunoff%flow(iunit,nt) = TRunoff%flow(iunit,nt) + TRunoff%erout(iunit,nt)
             enddo

             ! first round of trapping, for those main channel reservoirs that both regulate flow and trap sediment
             if(Tres_para%Eff_trapping(iunit)>TINYVALUE) then            
                 call res_trapping(iunit,nt_nmud)
                 Tres%wres(iunit,nt_nmud) = Tres%wres(iunit,nt_nmud) + Tres%dwres(iunit,nt_nmud) * localDeltaT
                 call res_trapping(iunit,nt_nsan)
                 Tres%wres(iunit,nt_nsan) = Tres%wres(iunit,nt_nsan) + Tres%dwres(iunit,nt_nsan) * localDeltaT
             end if

             !! TODO: second round of trapping, for those main-channel reservoirs that trap sediment only
             !if(Tres_para%Eff_trapping_r(iunit)>TINYVALUE) then             
             !    call res_trapping_r(iunit,nt_nmud)
             !    Tres%wres(iunit,nt_nmud) = Tres%wres(iunit,nt_nmud) + Tres%dwres(iunit,nt_nmud) * localDeltaT
             !    call res_trapping_r(iunit,nt_nsan)
             !    Tres%wres(iunit,nt_nsan) = Tres%wres(iunit,nt_nsan) + Tres%dwres(iunit,nt_nsan) * localDeltaT
             !end if

             do nt=nt_nmud,nt_nsan
               TRunoff%erowm_regf(iunit,nt) = TRunoff%erowm_regf(iunit,nt) - TRunoff%erout(iunit,nt)
               TRunoff%flow(iunit,nt) = TRunoff%flow(iunit,nt) - TRunoff%erout(iunit,nt)
             enddo
          end if


!#endif

       end do ! iunit
       endif  ! euler_calc     

       if (inundflag) then
            ! Channel -- floodplain exchange computation :      
              call ChnlFPexchg ( )
            ! update variables after channel-floodplain exchanges
            ! Floodplain water volume :
              TRunoff%wf_ini = TRunoff%wf_exchg
            ! Floodplain max water depth :
              TRunoff%hf_ini = TRunoff%hf_exchg     
            ! Floodplain area fraction (not including channel) ! will be deleted
              TRunoff%ff_ini = TRunoff%ff_fp
            ! Flooded area fraction (including channel):
              TRunoff%ffunit_ini = TRunoff%ff_unit
            ! Channel water depth
              TRunoff%yr(:,1) = TRunoff%yr_exchg
            ! Channel storage
              TRunoff%wr(:,1) = TRunoff%wr_exchg
            ! Aggregate net floodplain storage change from subcycle to timestep 
              TRunoff%se_rf = TRunoff%se_rf + TRunoff%netchange
        end if                   
       negchan = min(negchan, minval(TRunoff%wr(:,:)))
       call t_stopf('mosartr_chanroute') 
    end do  ! DLevelH2R
    ! subcycling within MOSART ends

   ! check for negative channel storage
    if (negchan < -1.e-10 .and. negchan >= -1.e-8) then
       write(iulog,*) 'Warning: Small negative channel storage found! ',negchan
    elseif(negchan < -1.e-8) then
       write(iulog,*) 'Error: Negative channel storage found! ',negchan
       call shr_sys_abort('mosart: negative channel storage')
    endif
    TRunoff%flow = TRunoff%flow / Tctl%DLevelH2R
    TRunoff%erowm_regi(:,nt_nmud:nt_nsan) = TRunoff%erowm_regi(:,nt_nmud:nt_nsan) / Tctl%DLevelH2R
    TRunoff%erowm_regf(:,nt_nmud:nt_nsan) = TRunoff%erowm_regf(:,nt_nmud:nt_nsan) / Tctl%DLevelH2R
    TRunoff%eroup_lagi = TRunoff%eroup_lagi / Tctl%DLevelH2R
    TRunoff%eroup_lagf = TRunoff%eroup_lagf / Tctl%DLevelH2R
    TRunoff%eroutup_avg = TRunoff%eroutup_avg / Tctl%DLevelH2R
    TRunoff%erlat_avg = TRunoff%erlat_avg / Tctl%DLevelH2R
    if (heatflag) then
       THeat%Ha_eroutup_avg = THeat%Ha_eroutup_avg / Tctl%DLevelH2R
       THeat%Ha_erlat_avg = THeat%Ha_erlat_avg / Tctl%DLevelH2R
       THeat%Tt_avg = THeat%Tt_avg / Tctl%DLevelH2R
       THeat%Tr_avg = THeat%Tr_avg / Tctl%DLevelH2R
    end if
    TRunoff%etexch_avg = TRunoff%etexch_avg / Tctl%DLevelH2R
    TRunoff%erexch_avg = TRunoff%erexch_avg / Tctl%DLevelH2R

    !------------------
    ! WRM Regulation
    ! WRM ExtractionRegulatedFlow
    ! Do not update wr after regulation or extraction from reservoir release. 
    ! Because of the regulation, the wr might get to crazy uncontrolled values, 
    ! assume in this case wr is not changed. The storage in reservoir handles it.
    !------------------

    if (wrmflag) then
       if (ctlSubwWRM%RegulationFlag>0) then
          ! compute the erowm_reg terms and adjust the flow diagnostic
          do iunit=rtmCTL%begr,rtmCTL%endr
             TRunoff%erowm_regi(iunit,nt_nliq) = -TRunoff%erout(iunit,nt_nliq)
             TRunoff%flow(iunit,nt_nliq) = TRunoff%flow(iunit,nt_nliq) + TRunoff%erout(iunit,nt_nliq)
          enddo
          localDeltaT = Tctl%DeltaT
!          call t_startf('mosartr_wrm_Reg')
!          do iunit=rtmCTL%begr,rtmCTL%endr
!             if (TUnit%mask(iunit) > 0) then
!                call Regulation(iunit, localDeltaT) !move regulation back into the subcycling, Tian 9/26/2018
!             endif
!          enddo
!          call t_stopf('mosartr_wrm_Reg')
          if (ctlSubwWRM%ExtractionFlag > 0 ) then
             call t_startf('mosartr_wrm_ERFlow')
             call ExtractionRegulatedFlow(localDeltaT)
             ! a simple treatment after extracting water from the regulated streamflow. Assuming the extraction won't change the water temperature in the release
             ! but the heat flux will be changed due to changing streamflow
             if (heatflag) then
                 do iunit=rtmCTL%begr,rtmCTL%endr
                     THeat%Ha_rout(iunit) = -cr_advectheat(abs(TRunoff%erout(iunit,nt_nliq)), THeat%Tr(iunit))
                 enddo
             end if
             call t_stopf('mosartr_wrm_ERFlow')
          endif
!          !--- now subtract updated erout to update flow calc
!          ! compute the erowm_reg terms and adjust the flow diagnostic
          do iunit=rtmCTL%begr,rtmCTL%endr
             TRunoff%erowm_regf(iunit,nt_nliq) = -TRunoff%erout(iunit,nt_nliq)
             TRunoff%flow(iunit,nt_nliq) = TRunoff%flow(iunit,nt_nliq) - TRunoff%erout(iunit,nt_nliq)
          enddo
       endif
    endif

    !------------------
    ! WRM post Euler updates
    !------------------

    if (wrmflag) then
       call t_startf('mosartr_wrm_estrfdef')
       call estimate_returnflow_deficit()
       if (ctlSubwWRM%ExtractionFlag > 0) then
          StorWater%deficit = StorWater%demand
          StorWater%deficit = StorWater%deficit/Tctl%DeltaT !change it back from mm/timestep to mm/s for output
       endif
       call t_stopf('mosartr_wrm_estrfdef')
    endif

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

    if(TUnit%tlen(iunit) <= TUnit%hlen(iunit)) then ! if no tributaries, not subnetwork channel routing
        TRunoff%etout(iunit,nt) = -TRunoff%etin(iunit,nt)
    else
        if(nt == nt_nliq) then
    !   !     !TRunoff%vt(iunit,nt) = CRVRMAN(TUnit%tslp(iunit), TUnit%nt(iunit), TRunoff%rt(iunit,nt))
            TRunoff%vt(iunit,nt) = CRVRMAN_nosqrt(TUnit%tslpsqrt(iunit), TUnit%nt(iunit), TRunoff%rt(iunit,nt))
            TRunoff%etout(iunit,nt) = -TRunoff%vt(iunit,nt) * TRunoff%mt(iunit,nt)
            if(TRunoff%wt(iunit,nt) + (TRunoff%etin(iunit,nt) + TRunoff%etout(iunit,nt)) * theDeltaT < TINYVALUE) then
              TRunoff%etout(iunit,nt) = -(TRunoff%etin(iunit,nt) + TRunoff%wt(iunit,nt)/theDeltaT)
              if(TRunoff%mt(iunit,nt) > 0._r8) then
                 TRunoff%vt(iunit,nt) = -TRunoff%etout(iunit,nt)/TRunoff%mt(iunit,nt)
              end if
            end if
        else
            TRunoff%etout(iunit,nt) = TRunoff%conc_t(iunit,nt)*TRunoff%etout(iunit,nt_nliq)
            if(TRunoff%etout(iunit,nt) < -TINYVALUE .and. &
               TRunoff%wt(iunit,nt) + (TRunoff%etin(iunit,nt) + TRunoff%etout(iunit,nt)) * theDeltaT < TINYVALUE) then
              TRunoff%etout(iunit,nt) = -(TRunoff%etin(iunit,nt) + TRunoff%wt(iunit,nt)/theDeltaT)
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

    if(Tctl%RoutingMethod == KW) then
       call Routing_KW(iunit, nt, theDeltaT)
    else if(Tctl%RoutingMethod == DW) then
       if ( use_ocn_rof_two_way ) then
          call Routing_DW_ocn_rof_two_way(iunit, nt, theDeltaT)
       else
          call Routing_DW(iunit, nt, theDeltaT)
       end if
    else
       call shr_sys_abort('Wrong routing method! There are only 2 methods available. 1==KW, 2==DW.')
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

    TRunoff%erin(iunit,nt) = TRunoff%erin(iunit,nt) - TRunoff%eroutUp(iunit,nt)

    ! estimate the outflow
    if(TUnit%rlen(iunit) <= 0._r8) then ! no river network, no channel routing
       TRunoff%vr(iunit,nt) = 0._r8
       TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
    else
       ! skip the channel routing if possible numerical instability
       if(TUnit%areaTotal2(iunit)/TUnit%rwidth(iunit)/TUnit%rlen(iunit) > 1e6_r8) then
          TRunoff%vr(iunit,nt) = 0._r8
          TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
       else
          if(nt == nt_nliq) then
              !TRunoff%vr(iunit,nt) = CRVRMAN(TUnit%rslp(iunit), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
              TRunoff%vr(iunit,nt) = CRVRMAN_nosqrt(TUnit%rslpsqrt(iunit), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
              TRunoff%erout(iunit,nt) = -TRunoff%vr(iunit,nt) * TRunoff%mr(iunit,nt)
              if(-TRunoff%erout(iunit,nt) > TINYVALUE .and. TRunoff%wr(iunit,nt) + &
                 (TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT < TINYVALUE) then
                 TRunoff%erout(iunit,nt) = -(TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%wr(iunit,nt)*MaxStorageDepleted/ theDeltaT)
                 if(TRunoff%mr(iunit,nt) > 0._r8) then
                    TRunoff%vr(iunit,nt) = -TRunoff%erout(iunit,nt) / TRunoff%mr(iunit,nt)
                 end if
              end if
          else
              TRunoff%erout(iunit,nt) = TRunoff%conc_r(iunit,nt) * TRunoff%erout(iunit,nt_nliq)
              if(-TRunoff%erout(iunit,nt) > TINYVALUE .and. TRunoff%wr(iunit,nt) + &
                 (TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT < TINYVALUE) then
                 TRunoff%erout(iunit,nt) = -(TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%wr(iunit,nt)*MaxStorageDepleted/ theDeltaT)
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

    integer  :: k, myflag 
    real(r8) :: temp_gwl, temp_dwr, temp_gwl0
    real(r8) :: temp1, temp2, w_temp
    real(r8 ) :: y_c, len_c, slp_c             ! Water depth (m), length (m) and bed slope (dimensionless) of the current channel.
    real(r8 ) :: y_down, len_down, slp_down    ! Water depth (m), length (m) and bed slope (dimensionless) of the downstream channel.
    
    ! estimate the inflow from upstream units
    TRunoff%erin(iunit,nt) = 0._r8

    TRunoff%erin(iunit,nt) = TRunoff%erin(iunit,nt) - TRunoff%eroutUp(iunit,nt)
    !do k=1, rtmCTL%nUp(iunit)
    !    TRunoff%erin(iunit,nt) = TRunoff%erin(iunit,nt) - TRunoff%erout(rtmCTL%iUp(iunit,k),nt)
    !enddo

    ! estimate the outflow
    if(TUnit%rlen(iunit) <= 0._r8) then ! no river network, no channel routing
       TRunoff%vr(iunit,nt) = 0._r8
       TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
    elseif(TUnit%areaTotal2(iunit)/TUnit%rwidth(iunit)/TUnit%rlen(iunit) > 1e6_r8) then ! skip the channel routing if possible numerical instability
       TRunoff%vr(iunit,nt) = 0._r8
       TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
    else
       !TODO. If this channel is at basin outlet (downstream is ocean), use the KW method
       if(rtmCTL%mask(iunit) .eq. 3) then 
          call Routing_KW(iunit, nt, theDeltaT)
       else
          if(nt == nt_nliq) then 

              if(TRunoff%rslp_energy(iunit) >= TINYVALUE) then ! flow is from current channel to downstream
                TRunoff%vr(iunit,nt) = CRVRMAN(TRunoff%rslp_energy(iunit), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
                TRunoff%erout(iunit,nt) = -TRunoff%vr(iunit,nt) * TRunoff%mr(iunit,nt)
                if(TRunoff%erin(iunit,nt)*theDeltaT + TRunoff%wr(iunit,nt) <= TINYVALUE) then! much negative inflow from upstream, 
                   TRunoff%vr(iunit,nt) = 0._r8
                   TRunoff%erout(iunit,nt) = 0._r8
                elseif(TRunoff%erout(iunit,nt) <= -TINYVALUE .and. TRunoff%wr(iunit,nt) + &
                   (TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT < TINYVALUE) then
                   TRunoff%erout(iunit,nt) = -(TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%wr(iunit,nt)*MaxStorageDepleted / theDeltaT)
                   if(TRunoff%mr(iunit,nt) > TINYVALUE) then
                      TRunoff%vr(iunit,nt) = -TRunoff%erout(iunit,nt) / TRunoff%mr(iunit,nt)
                   end if
                end if
              elseif(TRunoff%rslp_energy(iunit) <= -TINYVALUE) then ! flow is from downstream to current channel
                 TRunoff%vr(iunit,nt) = -CRVRMAN(abs(TRunoff%rslp_energy(iunit)), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
                 TRunoff%erout(iunit,nt) = -TRunoff%vr(iunit,nt) * TRunoff%mr(iunit,nt)
                 if(rtmCTL%nUp_dstrm(iunit) > 1) then
                     if(TRunoff%erin_dstrm(iunit,nt)*theDeltaT + TRunoff%wr_dstrm(iunit,nt)/rtmCTL%nUp_dstrm(iunit) <= TINYVALUE) then! much negative inflow from upstream,
                         TRunoff%vr(iunit,nt) = 0._r8
                         TRunoff%erout(iunit,nt) = 0._r8
                     elseif(TRunoff%erout(iunit,nt) >= TINYVALUE .and. TRunoff%wr_dstrm(iunit,nt)/rtmCTL%nUp_dstrm(iunit)- TRunoff%erout(iunit,nt) * theDeltaT < TINYVALUE) then
                        TRunoff%erout(iunit,nt) = TRunoff%wr_dstrm(iunit,nt)*MaxStorageDepleted / theDeltaT / rtmCTL%nUp_dstrm(iunit)
                       if(TRunoff%mr(iunit,nt) > TINYVALUE) then
                           TRunoff%vr(iunit,nt) = -TRunoff%erout(iunit,nt) / TRunoff%mr(iunit,nt)
                        end if
                     end if    
                 else
                     if(TRunoff%erin_dstrm(iunit,nt)*theDeltaT + TRunoff%wr_dstrm(iunit,nt) <= TINYVALUE) then! much negative inflow from upstream,
                         TRunoff%vr(iunit,nt) = 0._r8
                         TRunoff%erout(iunit,nt) = 0._r8
                     elseif(TRunoff%erout(iunit,nt) >= TINYVALUE .and. TRunoff%wr_dstrm(iunit,nt) &
                       - TRunoff%erout(iunit,nt) * theDeltaT < TINYVALUE) then
                        TRunoff%erout(iunit,nt) = TRunoff%wr_dstrm(iunit,nt)*MaxStorageDepleted / theDeltaT
                       if(TRunoff%mr(iunit,nt) > TINYVALUE) then
                           TRunoff%vr(iunit,nt) = -TRunoff%erout(iunit,nt) / TRunoff%mr(iunit,nt)
                        end if
                     end if
                 end if                 
                 !TRunoff%vr(iunit,nt) = 0._r8
                 !TRunoff%erout(iunit,nt) = 0._r8
              else  ! no flow between current channel and downstream
                TRunoff%vr(iunit,nt) = 0._r8
                TRunoff%erout(iunit,nt) = 0._r8
              end if
          else
            if(TRunoff%erout(iunit,nt_nliq) <= -TINYVALUE) then ! flow is from current channel to downstream
              TRunoff%erout(iunit,nt) = TRunoff%conc_r(iunit,nt) * TRunoff%erout(iunit,nt_nliq)
              if(TRunoff%erin(iunit,nt)*theDeltaT + TRunoff%wr(iunit,nt) <= TINYVALUE) then! much negative inflow from upstream, 
                 TRunoff%erout(iunit,nt) = 0._r8
              elseif(TRunoff%erout(iunit,nt) <= -TINYVALUE .and. TRunoff%wr(iunit,nt) + &
                 (TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT < TINYVALUE) then
                 TRunoff%erout(iunit,nt) = -(TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%wr(iunit,nt)*MaxStorageDepleted / theDeltaT)
              end if

            elseif(TRunoff%erout(iunit,nt_nliq) >= TINYVALUE) then ! flow is from downstream to current channel
              TRunoff%erout(iunit,nt) = 0._r8
            else
              TRunoff%erout(iunit,nt) = 0._r8
            end if
          end if
       end if  
    end if
    
    if(TRunoff%erin(iunit,nt) < -TINYVALUE .and. TRunoff%erout(iunit,nt) < -TINYVALUE) then
        if((TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT + TRunoff%wr(iunit,nt) < 0._r8) then
            TRunoff%erout(iunit,nt) = 0._r8
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

    !if(iunit==103833) then
    !    write(unit=2110,fmt="(i10,9(e12.3))") myflag, TRunoff%wr(iunit,1), TRunoff%dwr(iunit,nt), TRunoff%erlateral(iunit,1), TRunoff%erin(iunit,1), TRunoff%erout(iunit,1), TRunoff%rslp_energy(iunit),TUnit%rslp(iunit), TRunoff%vr(iunit, nt_nliq), TRunoff%yr(iunit, nt_nliq)
    !     write(unit=4111,fmt="(i10, 7(e12.3))") myflag, TRunoff%rslp_energy(iunit), TRunoff%vr(iunit, nt_nliq),TUnit%rslp(iunit), TRunoff%erin(iunit,1), TRunoff%erout(iunit,1), TUnit%nr(iunit), TRunoff%rr(iunit,nt)
    !    write(unit=2112,fmt="(4(i10), 2(e12.3))") myflag,rtmCTL%mask(iunit), rtmCTL%nUp(iunit), rtmCTL%iUp(iunit,1), TUnit%rlen(iunit), TRunoff%erout(rtmCTL%iUp(iunit,1),1) 
    !end if

! check for stability
!    if(TRunoff%vr(iunit,nt) < -TINYVALUE .or. TRunoff%vr(iunit,nt) > 30) then
!       write(iulog,*) "Numerical error inRouting_DW, ", iunit,nt,TRunoff%vr(iunit,nt)
!    end if

 !check for negative wr
   ! if((TRunoff%wr(iunit,nt)/theDeltaT + TRunoff%dwr(iunit,nt))/TRunoff%wr(iunit,nt) < -TINYVALUE) then
   !    write(iulog,*) 'negative wr! -- Routing_DW', iunit, TRunoff%wr(iunit,nt), TRunoff%erlateral(iunit,nt), TRunoff%erin(iunit,nt), TRunoff%erout(iunit,nt), temp_gwl
       !stop          
   ! end if     
   
  end subroutine Routing_DW

!-----------------------------------------------------------------------

  subroutine Routing_DW_ocn_rof_two_way(iunit, nt, theDeltaT)
  ! !DESCRIPTION: diffusion wave routing method that considers coupling with ocn model
    implicit none
    integer, intent(in) :: iunit, nt
    real(r8), intent(in) :: theDeltaT
    character(len=*),parameter :: subname = '(Routing_DW_ocn_rof_two_way)'

    integer  :: k, myflag
    real(r8) :: temp_gwl, temp_dwr, temp_gwl0
    real(r8) :: temp1, temp2, w_temp
    real(r8 ) :: y_c, len_c, slp_c             ! Water depth (m), length (m) and bed slope (dimensionless) of the current channel.
    real(r8 ) :: y_down, len_down, slp_down    ! Water depth (m), length (m) and bed slope (dimensionless) of the downstream channel.

    ! estimate the inflow from upstream units
    TRunoff%erin(iunit,nt) = 0._r8

    TRunoff%erin(iunit,nt) = TRunoff%erin(iunit,nt) - TRunoff%eroutUp(iunit,nt)

    ! estimate the outflow
    if(TUnit%rlen(iunit) <= 0._r8) then ! no river network, no channel routing
       TRunoff%vr(iunit,nt) = 0._r8
       TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
    elseif(TUnit%areaTotal2(iunit)/TUnit%rwidth(iunit)/TUnit%rlen(iunit) > 1e6_r8) then ! skip the channel routing if possible numerical instability
       TRunoff%vr(iunit,nt) = 0._r8
       TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
    else
       ! when ocn rof two-way coupling is on, use DW in the specified outlets
       if ( .not. use_ocn_rof_two_way .and. rtmCTL%mask(iunit) .eq. 3 ) then !If this channel is at basin outlet (downstream is ocean), use the KW method
          call Routing_KW(iunit, nt, theDeltaT)
       elseif ( use_ocn_rof_two_way .and. rtmCTL%mask(iunit) .eq. 3 .and. TUnit%ocn_rof_coupling_ID(iunit) .eq. 0 ) then
          call Routing_KW(iunit, nt, theDeltaT)
       else
!          TODO: conc_r
!          if(nt == nt_nliq) then
              if(TRunoff%rslp_energy(iunit) >= TINYVALUE) then ! flow is from current channel to downstream
                TRunoff%vr(iunit,nt) = CRVRMAN(TRunoff%rslp_energy(iunit), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
                TRunoff%erout(iunit,nt) = -TRunoff%vr(iunit,nt) * TRunoff%mr(iunit,nt)
                if(TRunoff%erin(iunit,nt)*theDeltaT + TRunoff%wr(iunit,nt) <= TINYVALUE) then! much negative inflow from upstream, 
                   TRunoff%vr(iunit,nt) = 0._r8
                   TRunoff%erout(iunit,nt) = 0._r8
                elseif(TRunoff%erout(iunit,nt) <= -TINYVALUE .and. TRunoff%wr(iunit,nt) + &
                   (TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT < TINYVALUE) then
                   TRunoff%erout(iunit,nt) = -(TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%wr(iunit,nt)*0.95_r8 / theDeltaT)
                   if(TRunoff%mr(iunit,nt) > TINYVALUE) then
                      TRunoff%vr(iunit,nt) = -TRunoff%erout(iunit,nt) / TRunoff%mr(iunit,nt)
                   end if
                end if
              elseif(TRunoff%rslp_energy(iunit) <= -TINYVALUE) then ! flow is from downstream to current channel
                 TRunoff%vr(iunit,nt) = -CRVRMAN(abs(TRunoff%rslp_energy(iunit)), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
                 TRunoff%erout(iunit,nt) = -TRunoff%vr(iunit,nt) * TRunoff%mr(iunit,nt)
                 ! two-way coupling update: allow flow from ocn
                 if(rtmCTL%nUp_dstrm(iunit) > 1) then
                     if( TUnit%ocn_rof_coupling_ID(iunit) .ne. 1 .and. TRunoff%erin_dstrm(iunit,nt)*theDeltaT + TRunoff%wr_dstrm(iunit,nt)/rtmCTL%nUp_dstrm(iunit) <= TINYVALUE) then! much negative inflow from upstream,
                         TRunoff%vr(iunit,nt) = 0._r8
                         TRunoff%erout(iunit,nt) = 0._r8
                     elseif( TUnit%ocn_rof_coupling_ID(iunit) .ne. 1 .and. TRunoff%erout(iunit,nt) >= TINYVALUE .and. TRunoff%wr_dstrm(iunit,nt)/rtmCTL%nUp_dstrm(iunit)- TRunoff%erout(iunit,nt) * theDeltaT < TINYVALUE) then
                         TRunoff%erout(iunit,nt) = TRunoff%wr_dstrm(iunit,nt)*0.95_r8 / theDeltaT / rtmCTL%nUp_dstrm(iunit)
                         if(TRunoff%mr(iunit,nt) > TINYVALUE) then
                            TRunoff%vr(iunit,nt) = -TRunoff%erout(iunit,nt) / TRunoff%mr(iunit,nt)
                         end if
                     end if
                 else
                     if( TUnit%ocn_rof_coupling_ID(iunit) .ne. 1 .and. TRunoff%erin_dstrm(iunit,nt)*theDeltaT + TRunoff%wr_dstrm(iunit,nt) <= TINYVALUE) then! much negative inflow from upstream,
                        TRunoff%vr(iunit,nt) = 0._r8
                        TRunoff%erout(iunit,nt) = 0._r8
                     elseif( TUnit%ocn_rof_coupling_ID(iunit) .ne. 1 .and. TRunoff%erout(iunit,nt) >= TINYVALUE .and. TRunoff%wr_dstrm(iunit,nt) &
                       - TRunoff%erout(iunit,nt) * theDeltaT < TINYVALUE) then
                        TRunoff%erout(iunit,nt) = TRunoff%wr_dstrm(iunit,nt)*0.95_r8 / theDeltaT
                        if(TRunoff%mr(iunit,nt) > TINYVALUE) then
                           TRunoff%vr(iunit,nt) = -TRunoff%erout(iunit,nt) / TRunoff%mr(iunit,nt)
                        end if
                     end if
                 end if
              else  ! no flow between current channel and downstream
                TRunoff%vr(iunit,nt) = 0._r8
                TRunoff%erout(iunit,nt) = 0._r8
              end if
       end if
    end if

    if(TRunoff%erin(iunit,nt) < -TINYVALUE .and. TRunoff%erout(iunit,nt) < -TINYVALUE) then
       if((TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT + TRunoff%wr(iunit,nt) < 0._r8) then
           TRunoff%erout(iunit,nt) = 0._r8
       end if
    end if

    temp_dwr = TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)
    temp_gwl = TRunoff%qgwl(iunit,nt) * TUnit%area(iunit) * TUnit%frac(iunit)
    temp_gwl0 = temp_gwl
    if(abs(temp_gwl) <= TINYVALUE) then
       temp_gwl = 0._r8
    end if
    if(temp_gwl < -TINYVALUE) then
       write(iulog,*) 'mosart: ERROR temp_gwl negative',iunit,nt,TRunoff%qgwl(iunit,nt)
       call shr_sys_abort('mosart: ERROR temp_gwl negative ')
       if(TRunoff%wr(iunit,nt) < TINYVALUE) then
          temp_gwl = 0._r8
       else
          if(TRunoff%wr(iunit,nt)/theDeltaT + temp_dwr + temp_gwl < -TINYVALUE) then
             temp_gwl = -(temp_dwr + TRunoff%wr(iunit,nt) / theDeltaT)
          end if
       end if
    end if

    TRunoff%dwr(iunit,nt) = TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt) + temp_gwl

  end subroutine Routing_DW_ocn_rof_two_way

!-----------------------------------------------------------------------

  subroutine updateState_hillslope(iunit,nt)
  ! !DESCRIPTION: update the state variables at hillslope
    implicit none    
    integer, intent(in) :: iunit, nt
    character(len=*),parameter :: subname = '(updateState_hillslope)'

    if(nt==nt_nliq) then
        TRunoff%yh(iunit,nt) = TRunoff%wh(iunit,nt) !/ TUnit%area(iunit) / TUnit%frac(iunit) 
    else
        TRunoff%yh(iunit,nt) = 0._r8
    end if

  end subroutine updateState_hillslope

!-----------------------------------------------------------------------

  subroutine updateState_subnetwork(iunit,nt)
  ! !DESCRIPTION: update the state variables in subnetwork channel
    implicit none    
    integer, intent(in) :: iunit,nt
    character(len=*),parameter :: subname = '(updateState_subnetwork)'

    if(nt == nt_nliq) then
       if(TUnit%tlen(iunit) > 0._r8 .and. TRunoff%wt(iunit,nt) > 0._r8) then
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
    else
        if(TRunoff%wt(iunit,nt_nliq) >= TINYVALUE .and. TRunoff%wt(iunit,nt) >= TINYVALUE) then
            TRunoff%conc_t(iunit,nt) = TRunoff%wt(iunit,nt)/TRunoff%wt(iunit,nt_nliq)
        else
            TRunoff%conc_t(iunit,nt) = 0._r8
        end if
    end if
  end subroutine updateState_subnetwork

!-----------------------------------------------------------------------

  subroutine updateState_mainchannel(iunit, nt)
  ! !DESCRIPTION: update the state variables in main channel
    implicit none    
    integer, intent(in) :: iunit, nt
    character(len=*),parameter :: subname = '(updateState_mainchannel)'

    if(nt == nt_nliq) then
       if(TUnit%rlen(iunit) > 0._r8 .and. TRunoff%wr(iunit,nt) > 0._r8) then
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
    else   
        if(TRunoff%wr(iunit,nt_nliq) >= TINYVALUE .and. TRunoff%wr(iunit,nt) >= TINYVALUE) then
            TRunoff%conc_r(iunit,nt) = TRunoff%wr(iunit,nt)/TRunoff%wr(iunit,nt_nliq)
        else
            TRunoff%conc_r(iunit,nt) = 0._r8
        end if
    end if 
  end subroutine updateState_mainchannel

!-----------------------------------------------------------------------
    
  function CRRSLP(iunit_) result(rslp_)
  ! ! Function for calculating the water surface slope between the current and downstream grids
    implicit none
    integer, intent(in) :: iunit_ ! local index of current grid
    real(r8)             :: rslp_            ! rslp_ is  slope [-]
    
    real(r8 ) :: y_c, len_c, slp_c             ! Water depth (m), length (m) and bed slope (dimensionless) of the current channel.
    real(r8 ) :: y_down, len_down, slp_down    ! Water depth (m), length (m) and bed slope (dimensionless) of the downstream channel.
    character(len=*),parameter :: subname = '(CRRSLP)'

    y_c = TRunoff%yr(iunit_,nt_nliq)
    len_c = TUnit%rlen(iunit_)
    slp_c = TUnit%rslp(iunit_)
    y_down = TRunoff%yr_dstrm(iunit_)
    len_down = TUnit%rlen_dstrm(iunit_)
    slp_down = TUnit%rslp_dstrm(iunit_)
       
    ! Calculate water surface slope ( from current-channel surface mid-point to downstream-channel surface mid-point ) :
    rslp_ = (len_down * slp_down + len_c * slp_c + 2._r8 * y_c - 2._r8 * y_down) / (len_c + len_down)

  end function CRRSLP

!-----------------------------------------------------------------------
    
  function CRVRMAN(slp_, n_, rr_) result(v_)
  ! Function for calculating channel velocity according to Manning's equation.
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
  ! Function for calculating channel velocity according to Manning's equation.
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
  ! Function for overland from hillslope into the sub-network channels
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
  ! Function for estimate wetted channel area
    implicit none
    real(r8), intent(in) :: wr_, rlen_      ! storage of water, channel length
    real(r8)             :: mr_             ! wetted channel area
    character(len=*),parameter :: subname = '(GRMR)'
    
    mr_ = wr_ / rlen_

  end function GRMR
  
!-----------------------------------------------------------------------

  function GRHT(mt_, twid_) result(ht_)
  ! Function for estimating water depth assuming rectangular channel
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
  ! Function for estimating wetted perimeter assuming rectangular channel
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
  ! Function for estimating hydraulic radius
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
  ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
  ! here assuming the channel cross-section consists of three parts, from bottom to up,
  ! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
  ! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
  ! part 3 is a rectagular with the width rwid0
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
  ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
  ! here assuming the channel cross-section consists of three parts, from bottom to up,
  ! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
  ! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
  ! part 3 is a rectagular with the width rwid0
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

  subroutine Leapfrog(iunit_, nt_, deltaT_)
  ! !DESCRIPTION: Purpose: leapfrog method for channel routing
    implicit none
    integer, intent(in) :: iunit_, nt_    !
    real(r8), intent(in) :: deltaT_ ! 

    real(r8) :: wrtemp_, k1_, k2_
    real(r8) :: erout1_, erout2_
    
    
    wrtemp_ = TRunoff%wr(iunit_,nt_)
    call mainchannelRouting(iunit_,nt_,deltaT_)
    erout1_ = TRunoff%erout(iunit_,nt_)
    k1_ = TRunoff%dwr(iunit_,nt_)
    TRunoff%wr(iunit_,nt_) = TRunoff%wr(iunit_,nt_) + k1_ * deltaT_ * 0.5_r8
    call UpdateState_mainchannel(iunit_,nt_)
    call mainchannelRouting(iunit_,nt_,deltaT_)
    erout2_ = TRunoff%erout(iunit_,nt_)
    k2_ = TRunoff%dwr(iunit_,nt_)
    !TRunoff%dwr(iunit_,nt_) =k1_ * 0.75_r8 + k2_ * 0.25_r8
    TRunoff%wr(iunit_,nt_) = wrtemp_ + (k1_ * 0.75_r8 + k2_ * 0.25_r8) * deltaT_
    TRunoff%erout(iunit_,nt_) = erout1_ * 0.75_r8 + erout2_ * 0.25_r8
    call UpdateState_mainchannel(iunit_,nt_)
    
    !call mainchannelRouting(iunit_,nt_,deltaT_)

  end subroutine Leapfrog

!-----------------------------------------------------------------------


  subroutine Leapfrog_sed(iunit_, deltaT_)
  ! !DESCRIPTION: Purpose: leapfrog method for channel routing
    implicit none
    integer, intent(in) :: iunit_    !
    real(r8), intent(in) :: deltaT_ ! 

    real(r8) :: wrtemp_san, wrtemp_mud,k1_san,k1_mud,k2_san,k2_mud
    real(r8) :: erout1_mud, erout2_mud
    real(r8) :: erout1_san, erout2_san
    real(r8) :: wr_al_temp_san, wr_al_temp_mud,k1_al_san,k1_al_mud,k2_al_san,k2_al_mud
    real(r8) :: k1_ersb,k1_ermb,k2_ersb,k2_ermb

    wrtemp_mud = TRunoff%wr(iunit_,nt_nmud)
    wr_al_temp_mud = TRunoff%wr_al(iunit_,nt_nmud)
    wrtemp_san = TRunoff%wr(iunit_,nt_nsan)
    wr_al_temp_san = TRunoff%wr_al(iunit_,nt_nsan)

    call mainchannelSediment(iunit_,deltaT_)
    erout1_mud = TRunoff%erout(iunit_,nt_nmud)
    k1_ermb = TSedi%ermb_r(iunit_)
    k1_mud = TRunoff%dwr(iunit_,nt_nmud)
    TRunoff%wr(iunit_,nt_nmud) = TRunoff%wr(iunit_,nt_nmud) + k1_mud * deltaT_ * 0.5_r8
    call UpdateState_mainchannel(iunit_,nt_nmud)        
    erout1_san = TRunoff%erout(iunit_,nt_nsan)
    k1_ersb = TSedi%ersb_r(iunit_)
    k1_san = TRunoff%dwr(iunit_,nt_nsan)
    TRunoff%wr(iunit_,nt_nsan) = TRunoff%wr(iunit_,nt_nsan) + k1_san * deltaT_ * 0.5_r8
    call UpdateState_mainchannel(iunit_,nt_nsan)

    k1_al_mud = TRunoff%dwr_al(iunit_,nt_nmud)
    TRunoff%wr_al(iunit_,nt_nmud) = TRunoff%wr_al(iunit_,nt_nmud) + k1_al_mud * deltaT_ * 0.5_r8
    k1_al_san = TRunoff%dwr_al(iunit_,nt_nsan)
    TRunoff%wr_al(iunit_,nt_nsan) = TRunoff%wr_al(iunit_,nt_nsan) + k1_al_san * deltaT_ * 0.5_r8


    call mainchannelSediment(iunit_,deltaT_)    
    erout2_mud = TRunoff%erout(iunit_,nt_nmud)
    k2_ermb = TSedi%ermb_r(iunit_)
    k2_mud = TRunoff%dwr(iunit_,nt_nmud)
    TRunoff%wr(iunit_,nt_nmud) = wrtemp_mud + (k1_mud * 0.75_r8 + k2_mud * 0.25_r8) * deltaT_
    TRunoff%erout(iunit_,nt_nmud) = erout1_mud * 0.75_r8 + erout2_mud * 0.25_r8
    call UpdateState_mainchannel(iunit_,nt_nmud)    
    erout2_san = TRunoff%erout(iunit_,nt_nsan)
    k2_ersb = TSedi%ersb_r(iunit_)
    k2_san = TRunoff%dwr(iunit_,nt_nsan)
    TRunoff%wr(iunit_,nt_nsan) = wrtemp_san + (k1_san * 0.75_r8 + k2_san * 0.25_r8) * deltaT_
    TRunoff%erout(iunit_,nt_nsan) = erout1_san * 0.75_r8 + erout2_san * 0.25_r8
    call UpdateState_mainchannel(iunit_,nt_nsan)

    TSedi%ersb_r(iunit_) = k1_ersb*0.75_r8 + k2_ersb*0.25_r8
    TSedi%ermb_r(iunit_) = k1_ermb*0.75_r8 + k2_ermb*0.25_r8

    k2_al_mud = TRunoff%dwr_al(iunit_,nt_nmud)
    TRunoff%wr_al(iunit_,nt_nmud) = wr_al_temp_mud + (k1_al_mud * 0.75_r8 + k2_al_mud * 0.25_r8) * deltaT_
    k2_al_san = TRunoff%dwr_al(iunit_,nt_nsan)
    TRunoff%wr_al(iunit_,nt_nsan) = wr_al_temp_san + (k1_al_san * 0.75_r8 + k2_al_san * 0.25_r8) * deltaT_


    !call mainchannelSediment(iunit_,deltaT_)

  end subroutine Leapfrog_sed
  
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
  
  end subroutine printTest

!-----------------------------------------------------------------------

end MODULE MOSART_physics_mod

