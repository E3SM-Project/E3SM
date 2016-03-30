!
MODULE WRM_modules
! Description: core code of the WRM. 
! 
! Developed by Nathalie Voisin Feb 2012
! REVISION HISTORY:
!-----------------------------------------------------------------------

! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8, SHR_KIND_CL, SHR_KIND_IN
  use shr_const_mod  , only : SHR_CONST_REARTH, SHR_CONST_PI
  use RunoffMod      , only : Trunoff, Tctl, Tunit, rtmCTL
  use RtmVar         , only : iulog, barrier_timers
  use RtmSpmd        , only : iam, mpicom_rof, mastertask, masterproc, &
                              MPI_REAL8,MPI_INTEGER,MPI_CHARACTER,MPI_LOGICAL,MPI_MAX
  use RtmTimeManager , only : get_curr_date
#if (1 == 0)
  use MOSART_physics_mod, only : updatestate_hillslope,	updatestate_subnetwork,	&
                                 updatestate_mainchannel, hillsloperouting, &
                                 subnetworkrouting, mainchannelrouting
  use rof_cpl_indices, only : nt_rtm
  use WRM_returnflow , only : insert_returnflow_channel
#endif
  use WRM_type_mod   , only : ctlSubwWRM, WRMUnit, StorWater, &
                              aVect_wg, aVect_wd, sMatP_g2d, sMatP_d2g, &
                              gsMap_wg, gsMap_wd
  use rof_cpl_indices, only : nt_nliq
  use mct_mod
  use perf_mod       , only : t_startf, t_stopf
     
  implicit none
  private

  real(r8), parameter :: MYTINYVALUE = 1.0e-50_r8  ! double precisio variable has a significance of about 16 decimal digits

! !PUBLIC MEMBER FUNCTIONS:
#if (1 == 0)
  public Euler_WRM  
#endif
  public irrigationExtraction
  public irrigationExtractionSubNetwork
  public irrigationExtractionMainChannel
  public RegulationRelease
  public Regulation
  public ExtractionRegulatedFlow
  public WRM_storage_targets
                                           
!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------

#if (1 == 0)
  subroutine Euler_WRM
     ! !DESCRIPTION: solve the ODEs with Euler algorithm
     implicit none

     integer :: iunit, m, k   !local index
     real(r8) :: temp_erout(nt_rtm), localDeltaT
     character(len=*),parameter :: subname='(Euler_WRM)'

     if (Tctl%RoutingFlag == 1) then
        !print*, "Running WRM Euler"
        do iunit=rtmCTL%begr,rtmCTL%endr
           call hillslopeRouting(iunit, Tctl%DeltaT)
           Trunoff%wh(iunit,:) = Trunoff%wh(iunit,:) + Trunoff%dwh(iunit,:) * Tctl%DeltaT
           call UpdateState_hillslope(iunit)
           Trunoff%etin(iunit,:) = (-Trunoff%ehout(iunit,:) + Trunoff%qsub(iunit,:)) * TUnit%area(iunit) * TUnit%frac(iunit)
        end do
        ! exctraction from available surface runoff within the subw
        if (ctlSubwWRM%ExtractionFlag > 0) then
           call irrigationExtractionSubNetwork
        end if

        Trunoff%flow = 0._r8
        do m=1,Tctl%DLevelH2R
           do iunit=rtmCTL%begr,rtmCTL%endr
              Trunoff%erlateral(iunit,:) = 0._r8
              do k=1,TUnit%numDT_t(iunit)
                 localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_t(iunit)
                 call subnetworkRouting(iunit,localDeltaT)
                 Trunoff%wt(iunit,:) = Trunoff%wt(iunit,:) + Trunoff%dwt(iunit,:) * localDeltaT
                 call UpdateState_subnetwork(iunit)
                 Trunoff%erlateral(iunit,:) = Trunoff%erlateral(iunit,:)-Trunoff%etout(iunit,:)
              end do
              Trunoff%erlateral(iunit,:) = Trunoff%erlateral(iunit,:) / TUnit%numDT_t(iunit)
           end do

           do iunit=rtmCTL%begr,rtmCTL%endr
              if (TUnit%mask(iunit) > 0) then
                 temp_erout = 0._r8
                 do k=1,TUnit%numDT_r(iunit)
                    localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_r(iunit)
                    call mainchannelRouting(iunit,localDeltaT)
                    Trunoff%wr(iunit,:) = Trunoff%wr(iunit,:) + Trunoff%dwr(iunit,:) * localDeltaT
                    call UpdateState_mainchannel(iunit)
                    temp_erout(:) = temp_erout(:) + Trunoff%erout(iunit,:) ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
                 end do
                 temp_erout(:) = temp_erout(:) / TUnit%numDT_r(iunit)
                 Trunoff%erout(iunit,:) = temp_erout(:)
                 localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
                 if (ctlSubwWRM%ExtractionMainChannelFlag > 0 .AND. ctlSubwWRM%ExtractionFlag > 0 ) then
                    call IrrigationExtractionMainChannel(iunit, localDeltaT )
                    if (ctlSubwWRM%TotalDemandFlag > 0 .AND. ctlSubwWRM%ReturnFlowFlag > 0 ) then
                       call insert_returnflow_channel(iunit, localDeltaT )
                    endif
                    ! update main channel storage as well
                    temp_erout(:) = temp_erout(:) - Trunoff%erout(iunit,:) ! change in erout after regulation and extraction
                    Trunoff%dwr(iunit,:) =  temp_erout(:)
                    Trunoff%wr(iunit,:) = Trunoff%wr(iunit,:) + Trunoff%dwr(iunit,:) * localDeltaT
                    call UpdateState_mainchannel(iunit)
                 endif
                 if ( ctlSubwWRM%RegulationFlag>0 .and. WRMUnit%INVicell(iunit) > 0 .and. WRMUnit%MeanMthFlow(iunit,13) > 0.01_r8 ) then
                    call Regulation(iunit, localDeltaT)
                    if ( ctlSubwWRM%ExtractionFlag > 0 ) then
                       call ExtractionRegulatedFlow(iunit, localDeltaT)
                    endif
                 endif

                 Trunoff%flow(iunit,:) = Trunoff%flow(iunit,:) - Trunoff%erout(iunit,:)
                 ! do not update wr after regulation or extraction from reservoir release. Because of the regulation, the wr might get to crazy uncontrolled values, assume in this case wr is not changed. The storage in reservoir handles it.

              end if
           end do
           !print*, "flow ", m ,Trunoff%flow(137,:), Trunoff%erlateral(137,:), StorWater%supply(137)
        end do

     else   ! routing flag
        Trunoff%flow = 0._r8
        do m=1,Tctl%DLevelH2R
           do iunit=rtmCTL%begr,rtmCTL%endr
              Trunoff%erlateral(iunit,:) = (Trunoff%qsur(iunit,:) + Trunoff%qsub(iunit,:)) * TUnit%area(iunit) * TUnit%frac(iunit)
           end do
           ! exctraction from available surface runoff within the subw
           if (ctlSubwWRM%ExtractionFlag > 0) then
              call irrigationExtraction
           end if

           ! regulation  performed prior to routing in order to route the release properly
           !if ( ctlSubwWRM%RegulationFlag > 0 ) then
           !   call RegulationRelease
           !endif

           do iunit=rtmCTL%begr,rtmCTL%endr
              temp_erout = 0._r8
              if (ctlSubwWRM%ExtractionMainChannelFlag > 0 .AND. ctlSubwWRM%ExtractionFlag > 0 ) then
                 localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
                 call IrrigationExtractionMainChannel(iunit, localDeltaT )
                 !print*,"done with extraction - up[date mainchannel ",iunit
                 if (ctlSubwWRM%TotalDemandFlag > 0 .AND. ctlSubwWRM%ReturnFlowFlag > 0 ) then
                    call insert_returnflow_channel(iunit, localDeltaT )
                 endif
                 call UpdateState_mainchannel(iunit)
              endif

              do k=1,TUnit%numDT_r(iunit)
                 localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_r(iunit)
                 call mainchannelRouting(iunit,localDeltaT)
                 Trunoff%wr(iunit,:) = Trunoff%wr(iunit,:) + Trunoff%dwr(iunit,:) * localDeltaT
                 call UpdateState_mainchannel(iunit)
                 temp_erout(:) = temp_erout(:) + Trunoff%erout(iunit,:) ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
              end do

              temp_erout(:) = temp_erout(:) / TUnit%numDT_r(iunit)
              Trunoff%erout(iunit,:) = temp_erout(:)
              localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
              !if (ctlSubwWRM%ExtractionMainChannelFlag > 0 .AND. ctlSubwWRM%ExtractionFlag > 0 ) then
              !   call IrrigationExtractionMainChannel(iunit, localDeltaT )
              !   if (ctlSubwWRM%TotalDemandFlag > 0 .AND. ctlSubwWRM%ReturnFlowFlag > 0 ) then
              !     call insert_returnflow_channel(iunit, localDeltaT )
              !   endif
              !   ! update main channel storage as well
              !   temp_erout(:) = temp_erout(:) - Trunoff%erout(iunit,:) ! change in erout after regulation and extraction
              !   Trunoff%dwr(iunit,:) =  temp_erout(:)
              !   Trunoff%wr(iunit,:) = Trunoff%wr(iunit,:) + Trunoff%dwr(iunit,:) * localDeltaT
              !   call UpdateState_mainchannel(iunit)
              !endif

              if ( ctlSubwWRM%RegulationFlag>0 .and. WRMUnit%INVicell(iunit) > 0 ) then
                 call Regulation(iunit, localDeltaT)
                 if ( ctlSubwWRM%ExtractionFlag > 0 .and. WRMUnit%subw_Ndepend(iunit) > 0) then
                    call ExtractionRegulatedFlow(iunit, localDeltaT)
                 endif
              endif
              Trunoff%flow(iunit,:) = Trunoff%flow(iunit,:) - Trunoff%erout(iunit,:)

           end do
        end do
     end if

     Trunoff%flow = Trunoff%flow / Tctl%DLevelH2R

  end subroutine Euler_WRM
#endif
!-----------------------------------------------------------------------

  subroutine irrigationExtraction
     ! !DESCRIPTION: subnetwork channel routing irrigation extraction
     implicit none    
     integer :: iunit      ! local index
     real(r8) :: flow_vol, temp         ! flow in cubic meter rather than cms
     character(len=*),parameter :: subname='(irrigationExtraction)'

     do iunit=rtmCTL%begr,rtmCTL%endr

        flow_vol = Trunoff%erlateral(iunit,nt_nliq) * Tctl%DeltaT/Tctl%DLevelH2R
        temp = flow_vol
        if (TUnit%mask(iunit) > 0) then
           if ( flow_vol >= StorWater%demand(iunit) ) then 
              StorWater%supply(iunit)= StorWater%supply(iunit) + StorWater%demand(iunit)
              flow_vol = flow_vol - StorWater%demand(iunit)
              StorWater%demand(iunit)= 0._r8
           else
              StorWater%supply(iunit)= StorWater%supply(iunit) + flow_vol
              StorWater%demand(iunit)= StorWater%demand(iunit) - flow_vol
              flow_vol = 0._r8
           end if 
! dwt is not updated because extraction is taken from water getting out of subnetwork channel routing only
        end if
        Trunoff%erlateral(iunit,nt_nliq) = flow_vol / (Tctl%DeltaT/Tctl%DLevelH2R)
        !if ( StorWater%demand(iunit) > 0 .and. iunit.eq.96) then
        !  print*, temp, flow_vol, StorWater%supply(iunit), StorWater%demand(iunit)
        !endif 
     end do
  end subroutine irrigationExtraction

!-----------------------------------------------------------------------

  subroutine irrigationExtractionSubNetwork
     ! !DESCRIPTION: subnetwork channel routing irrigation extraction
     implicit none
     integer :: iunit      ! local index
     real(r8) :: flow_vol, temp, temp_vol         ! flow in cubic meter rather than cms
     character(len=*),parameter :: subname='(irrigationExtractionSubNetwork)'

     do iunit=rtmCTL%begr,rtmCTL%endr

        flow_vol = Trunoff%etin(iunit,nt_nliq)* Tctl%DeltaT/Tctl%DLevelH2R
        if (TUnit%mask(iunit) > 0) then
           if ( flow_vol >= StorWater%demand(iunit) ) then
              StorWater%supply(iunit)= StorWater%supply(iunit) + StorWater%demand(iunit)
              flow_vol = flow_vol - StorWater%demand(iunit)
              StorWater%demand(iunit)= 0._r8
           else
              StorWater%supply(iunit)= StorWater%supply(iunit) + flow_vol
              StorWater%demand(iunit)= StorWater%demand(iunit) - flow_vol
              flow_vol = 0._r8
           end if
        end if
        Trunoff%etin(iunit,nt_nliq) = flow_vol / (Tctl%DeltaT/Tctl%DLevelH2R)
     end do
  end subroutine irrigationExtractionSubNetwork

!-----------------------------------------------------------------------

  subroutine irrigationExtractionMainChannel(iunit, TheDeltaT )
     ! !DESCRIPTION: main channel routing irrigation extraction - restrict to 50% of the flow, something needs to flow else instability
     implicit none
     integer :: match
     integer, intent(in) :: iunit
     real(r8), intent(in) :: theDeltaT
     real(r8) :: flow_vol, frac         ! flow in cubic meter rather than cms
     character(len=*),parameter :: subname='(irrigationExtractionMainChannel)'

     match = 0
     frac = 0.5_r8 ! control the fraction of the flow that can be extracted

     ! added if statement for test
     if (Trunoff%wr(iunit,nt_nliq) > MYTINYVALUE .and. StorWater%demand(iunit) > MYTINYVALUE) then
        if (TUnit%mask(iunit) > 0 .and. TUnit%rlen(iunit) > MYTINYVALUE) then
           !flow_vol = -Trunoff%erout(iunit,nt_nliq) * theDeltaT
           flow_vol = Trunoff%wr(iunit,nt_nliq)

           if ( (frac*flow_vol) >= StorWater%demand(iunit) ) then
              StorWater%supply(iunit)= StorWater%supply(iunit) + StorWater%demand(iunit)
              flow_vol = flow_vol - StorWater%demand(iunit)
              StorWater%demand(iunit)= 0._r8
           else
              StorWater%supply(iunit)= StorWater%supply(iunit) + frac*flow_vol
              StorWater%demand(iunit)= StorWater%demand(iunit) - frac*flow_vol
              flow_vol = (1._r8-frac)*flow_vol
           end if
           !Trunoff%erout(iunit,nt_nliq) = -flow_vol / (theDeltaT)
           Trunoff%wr(iunit,nt_nliq) = flow_vol
           if ( Trunoff%wr(iunit,nt_nliq) < MYTINYVALUE ) then
              print*, "error with extraction from main chanel, ",iunit, Trunoff%wr(iunit,nt_nliq)
           endif
        endif
     endif

     !if (  match > 0 ) then
     !   print*, "MAIN after extract", StorWater%demand(137), StorWater%supply(137), -Trunoff%erout(137,nt_nliq)
     !endif

  end subroutine irrigationExtractionMainChannel

!-----------------------------------------------------------------------
  
  subroutine RegulationRelease
     !! DESCRIPTION: computes the expected monthly release based on Biemans (2011)
     implicit none
     integer :: idam, yr, month, day, tod
     real(r8) :: factor, k
     character(len=*),parameter :: subname='(RegulationRelease)'

     k = 1.0_r8

     call get_curr_date(yr, month, day, tod)

     do idam=1,ctlSubwWRM%LocalNumDam

        !if ( WRMUnit%use_FCon(idam) > 0 .or. WRMUnit%use_Supp(idam) > 0) then
           StorWater%release(idam) = WRMUnit%MeanMthFlow(idam,13)
        !endif
        ! prioity to irrigation and other use since storage targets for FC
        !if ( WRMUnit%use_Elec(idam) > 0 .or. WRMUnit%use_Irrig(idam) >0) then
           k = 1._r8
           factor = 0._r8
           k = WRMUnit%StorMthStOp(idam) / ( 0.85 * WRMUnit%StorCap(idam) )
           if ( WRMUnit%INVc(idam) .gt. 0.1_r8 ) then
              factor = (1._r8/(0.5_r8*WRMUnit%INVc(idam)))*(1._r8/(0.5_r8*WRMUnit%INVc(idam)))
           endif
              
           if ( WRMUnit%use_Elec(idam) > 0 .or. WRMUnit%use_Irrig(idam) >0) then
              !if ( (1._r8/WRMUnit%INVc(idam)) >= 0.5_r8 ) then
              if ( WRMUnit%INVc(idam) <= 2._r8 ) then
                 StorWater%release(idam) = k * Storwater%pre_release(idam,month)
              else
                 StorWater%release(idam) = k * factor*Storwater%pre_release(idam,month) + (1._r8-factor) * WRMUnit%MeanMthFlow(idam,month)  
              end if
           else
              if ( WRMUnit%INVc(idam) <= 2._r8 ) then
                 StorWater%release(idam) = k * WRMUnit%MeanMthFlow(idam, 13)
              else
                 StorWater%release(idam) = k * factor*WRMUnit%MeanMthFlow(idam, 13) + (1._r8-factor) * WRMUnit%MeanMthFlow(idam, month)
              endif
           !else
           !   StorWater%release(idam) = WRMUnit%MeanMthFlow(idam,month)
           !endif
           ! Run-on-the-river flow
           !if ( WRMUnit%use_FCon(idam) .eq. 0 .and.  WRMUnit%use_Irrig(idam).eq.0 .and. WRMUnit%use_Elec(idam) > 0 ) then
           !  StorWater%release(idam) = WRMUnit%MeanMthFlow(idam,month)
           end if
! PRIORITY TO INTEGRATION
!          if ( WRMUnit%use_FCon(idam) > 0 ) then
!             StorWater%release(idam) = WRMUnit%MeanMthFlow(idam,13)
!          endif

!       endif

     end do

  end subroutine RegulationRelease

!-----------------------------------------------------------------------

  subroutine WRM_storage_targets
     ! !DESCRIPTION: definr the necessary drop in storage based in sotrage at srta of the month
     ! NOT TO BE RUN IN EULER
     implicit none
     integer :: idam, yr, month, day, tod
     integer nio,ierror                 ! unit number of a file, flag number of IO status
     integer :: mth, Nmth, Nmth_fill       ! number of sign change
     real(r8) :: drop, fill, diff
     character(len=*),parameter :: subname='(WRM_storage_targets)'

     call get_curr_date(yr, month, day, tod)

     do idam=1,ctlSubwWRM%LocalNumDam

        drop = 0
        Nmth = 0
        if (WRMUnit%StorageCalibFlag(idam).eq.0) then
           if ( WRMUnit%use_FCon(idam) > 0 .and. WRMUnit%MthStFC(idam) > 0) then ! in the context of FC has priority
              ! modify release in order to maintain a certain storage level
              if ( WRMUnit%MthStFC(idam) <= WRMUnit%MthNdFC(idam) ) then
                 do mth = 1,12
                    if ( mth >= WRMUnit%MthStFC(idam) .and. mth < WRMUnit%MthNdFC(idam)) then
                       if ( WRMUnit%MeanMthFlow(idam, mth) >= WRMUnit%MeanMthFlow(idam, 13) ) then
                          drop = drop + 0._r8
                       else
                          drop = drop + abs(WRMUnit%MeanMthFlow(idam, 13) - WRMUnit%MeanMthFlow(idam, mth))
                       endif
                       Nmth = Nmth + 1
                    endif
                 enddo
              else if ( WRMUnit%MthStFC(idam) > WRMUnit%MthNdFC(idam) ) then
                 do mth =1,12
                    if (mth >= WRMUnit%MthStFC(idam) .or. mth < WRMUnit%MthNdFC(idam)) then
                       if ( WRMUnit%MeanMthFlow(idam, mth) >= WRMUnit%MeanMthFlow(idam, 13) ) then
                          drop = drop + 0._r8
                       else
                          drop = drop + abs(WRMUnit%MeanMthFlow(idam, 13) - WRMUnit%MeanMthFlow(idam, mth))
                       endif
                       Nmth = Nmth + 1
                    endif
                 enddo
              endif

              if ( Nmth > 0 ) then
                 if ( WRMUnit%MthStFC(idam) <= WRMUnit%MthNdFC(idam) ) then
                    if ( month >= WRMUnit%MthStFC(idam) .and. month < WRMUnit%MthNdFC(idam)) then
                       StorWater%release(idam) = StorWater%release(idam) + drop/Nmth
                    endif
                 else if ( WRMUnit%MthStFC(idam) > WRMUnit%MthNdFC(idam) ) then
                    if ( month >= WRMUnit%MthStFC(idam) .or. month < WRMUnit%MthNdFC(idam)) then
                       StorWater%release(idam) = StorWater%release(idam) + drop/Nmth
                    endif
                 endif
              endif

              ! now need to make sure that it will fill up but issue with spilling  in certain hydro-climatic conditions
              fill = 0  
              Nmth_fill = 0
              if ( WRMUnit%MthNdFC(idam) <= WRMUnit%MthStOP(idam) ) then
                 if ( month >= WRMUnit%MthNdFC(idam) .and. month < WRMUnit%MthStOp(idam) ) then
                    do mth = WRMUnit%MthNdFC(idam), WRMUnit%MthStOP(idam)
                       if ( WRMUnit%MeanMthFlow(idam, mth) > WRMUnit%MeanMthFlow(idam, 13) ) then
                          fill = fill + abs(WRMUnit%MeanMthFlow(idam, 13) - WRMUnit%MeanMthFlow(idam, mth))
                          Nmth_fill = Nmth_fill + 1
                       endif
                    end do
                    ! does drop fill up the reservoir?
                    !if ( fill > drop .and. Nmth_fill > 0 ) then
                    !   StorWater%release(idam) = WRMUnit%MeanMthFlow(idam, 13) + (fill - drop) / Nmth_fill
                    !else  !need to fill this reservoir 
                    if ( StorWater%release(idam) > WRMUnit%MeanMthFlow(idam, 13) ) then
                       StorWater%release(idam) = WRMUnit%MeanMthFlow(idam, 13)
                    endif
                    !endif
                 end if
              else if ( WRMUnit%MthNdFC(idam) > WRMUnit%MthStOP(idam) ) then
                 if ( month >= WRMUnit%MthNdFC(idam) .or. month < WRMUnit%MthStOp(idam)) then
                    do mth = WRMUnit%MthNdFC(idam), 12
                       if ( WRMUnit%MeanMthFlow(idam, mth) > WRMUnit%MeanMthFlow(idam, 13) ) then
                          fill = fill + abs(WRMUnit%MeanMthFlow(idam, 13) - WRMUnit%MeanMthFlow(idam, mth))
                          Nmth_fill = Nmth_fill + 1
                       endif
                    end do
                    do mth = 1, WRMUnit%MthStOP(idam)
                       if ( WRMUnit%MeanMthFlow(idam, mth) > WRMUnit%MeanMthFlow(idam, 13) ) then
                          fill = fill + abs(WRMUnit%MeanMthFlow(idam, 13) - WRMUnit%MeanMthFlow(idam, mth))
                          Nmth_fill = Nmth_fill + 1
                       endif
                    end do
                    ! does drop fill up the reservoir?
                    !if ( fill > drop .and. Nmth_fill > 0 ) then
                    !   StorWater%release(idam) = WRMUnit%MeanMthFlow(idam, 13) + (fill - drop) / Nmth_fill
                    !else  !need to fill this reservoir
                    if ( StorWater%release(idam) > WRMUnit%MeanMthFlow(idam, 13) ) then
                       StorWater%release(idam) = WRMUnit%MeanMthFlow(idam, 13)
                    endif
                    !endif
                 end if
              endif
           endif

           !!additional constraint when both FC and irrigation, but targets not computer. Inn previous version, priority given to irrigation. It wiorks great of the US but not when flooding is in su,m,mer - does not allow for multiple season growth
           !if ( WRMUnit%use_FCon(idam) > 0 .and. WRMUnit%MthStFC(idam) .eq. 0) then
  !if ( WRMUnit%use_Elec(idam) > 0 .or. WRMUnit%use_Irrig(idam) >0) then 
           !   StorWater%release(idam) = WRMUnit%MeanMthFlow(idam, 13)
           !endif
        else !enforce the strage targets from altimetry
           !calib4 code
           !fill = WRMUnit%MeanMthFlow(idam, month) + (StorWater%storage(idam)-WRMUnit%StorTarget(idam,month))/86400/30._r8
           mth = month+1
           if (mth>12) then
              mth = 1
           endif
           fill = (StorWater%storage(idam)-WRMUnit%StorTarget(idam,mth))/86400/30._r8
           !! for very large reservoir, keep releasing the annual mean monthly flow
           if ( WRMUnit%INVc(idam) <= 2._r8 ) then
              !fill = WRMUnit%MeanMthFlow(idam, 13) + (WRMUnit%StorTarget(idam,13)-WRMUnit%StorTarget(idam,month))/86400/30._r8
              !fill = StorWater%release(idam) + (WRMUnit%StorTarget(idam,13)-WRMUnit%StorTarget(idam,month))/86400/30._r8
              !fill = (WRMUnit%StorTarget(idam,13)-WRMUnit%StorTarget(idam,month))/86400/30._r8
              fill = (StorWater%storage(idam)-WRMUnit%StorTarget(idam,mth))/86400/30._r8
           endif
           !calib5 code
           !diff = StorWater%storage(idam) - WRMUnit%StorTarget(idam,month)
           !fill = WRMUnit%MeanMthFlow(idam, month) + (2*diff)/86400/30._r8
           ! for very large reservoir, keep releasing the annual mean monthly flow
           !if ( WRMUnit%INVc(idam) <= 2._r8 ) then
           !  fill = WRMUnit%MeanMthFlow(idam, 13) + (2*diff)/86400/30._r8
           !endif

           ! end change in code
           ! add constraints on releases 
           StorWater%release(idam) = fill 
!minimal constraint
           drop = 0.1_r8 * WRMUnit%MeanMthFlow(idam, month)
           StorWater%release(idam) = max(drop,fill)
!constraint option
           !if ( WRMUnit%INVc(idam) <= 2._r8 ) then
           !  drop = 0.5_r8*WRMUnit%MeanMthFlow(idam, 13)
           !else
           !  drop = 0.2_r8*WRMUnit%MeanMthFlow(idam, month)
           !endif   
           !  if ( fill .lt. drop ) then
           !     StorWater%release(idam) = drop
           !  endif

           !if ( WRMUnit%INVc(idam) <= 2._r8 ) then
           !  drop = max(1.5_r8*WRMUnit%MeanMthFlow(idam,13), WRMUnit%MeanMthFlow(idam, month))
           !else
           !  drop = max(2._r8*WRMUnit%MeanMthFlow(idam,13), WRMUnit%MeanMthFlow(idam, month))
           !endif
           !if ( fill .gt. drop) then 
           !      StorWater%release(idam) = drop
           !endif
           !print*,"calibrate release ",idam, month,fill, drop, WRMUnit%StorCap(idam)
           !print*,WRMUnit%MeanMthFlow(idam,month),WRMUnit%MeanMthFlow(idam,13),StorWater%storage(idam),WRMUnit%StorTarget(idam,month)
        endif
     end do

  end subroutine WRM_storage_targets

!-----------------------------------------------------------------------

  subroutine Regulation(iunit, TheDeltaT)
     !! DESCRIPTION: regulation of the flow from the reservoirs. The Regulation is applied to the flow entering the grid cell, i.e. the subw downstream of the reservoir. 
     ! !DESCRIPTION: CHANGE IN PLANS seems like erin get overwritten now play with erout
     implicit none
     integer , intent(in) :: iunit
     real(r8), intent(in) :: TheDeltaT
     !--- local ---
     integer  :: match
     integer  :: yr, month, day, tod
     integer  :: damID,k
     real(r8) :: flow_vol, flow_res, min_flow, min_stor, evap, max_stor
     character(len=*),parameter :: subname='(Regulation)'

     match = 0
     call get_curr_date(yr, month, day, tod)

     damID = WRMUnit%INVicell(iunit) 
     if ( damID > ctlSubwWRM%LocalNumDam .OR. damID <= 0 ) then
        print*, "Error in Regulation with DamID ",damID
        stop
     end if
     flow_vol = -Trunoff%erout(iunit,nt_nliq) * theDeltaT
     flow_res = StorWater%release(damID) * theDeltaT
     evap = StorWater%pot_evap(iunit) * theDeltaT * WRMUnit%SurfArea(damID) * 1000000._r8 ! potential evaporation in the grid cell the reservoir is
     min_flow = 0.1_r8 * WRMUnit%MeanMthFlow(damID, month) * theDeltaT
     min_stor = 0.1_r8 * WRMUnit%StorCap(damID)
     max_stor = WRMUnit%StorCap(damID)
     if ( WRMUnit%StorageCalibFlag(damID).eq.1) then
        ! code calib4
        min_stor = max(0.9_r8 * WRMUnit%MinStorTarget(damID),0.1_r8 * WRMUnit%StorCap(damID)) ! allows some variation
        ! code calib 5
        !k = WRMUnit%StorMthStOp(damID) / ( 0.85 * WRMUnit%StorCap(damID) )
        !min_stor = max(k * WRMUnit%MinStorTarget(damID),0.1_r8 * WRMUnit%StorCap(damID)) ! allows interannual variation
        ! end change in code 
        max_stor = min(1.1_r8*WRMUnit%MaxStorTarget(damID),WRMUnit%StorCap(damID))
     endif

     !print*, "preregulation", damID, StorWater%storage(damID), flow_vol, flow_res, min_flow, month

     if ( ( flow_vol + StorWater%storage(damID) - flow_res- evap) >= max_stor ) then
        flow_res =  flow_vol + StorWater%storage(damID) - max_stor - evap
        StorWater%storage(damID) = max_stor
     else if ( ( flow_vol + StorWater%storage(damID) - flow_res - evap) < min_stor ) then
        if ( flow_res<=(flow_vol-evap)) then
           StorWater%storage(damID) = StorWater%storage(damID) + flow_vol - flow_res - evap
        else if ( (flow_vol-evap)>=min_flow) then
           StorWater%storage(damID) = StorWater%storage(damID) 
           flow_res = flow_vol - evap
           !print*, "WARNING  No regulation", flow_vol, min_flow, damID
        else
           !print*, "ERROR - drying out the reservoir"
           !print*,damID, flow_vol, flow_res, StorWater%storage(damID), min_stor, min_flow
           !stop
           !flow_res = 0._r8
           !StorWater%storage(damID) = StorWater%storage(damID) - flow_res + flow_vol - evap
           flow_res = flow_vol
           StorWater%storage(damID) = StorWater%storage(damID) - flow_res + flow_vol - evap
           if (StorWater%storage(damID) < 0._r8) StorWater%storage(damID) = 0._r8
        endif
     else
        StorWater%storage(damID) = StorWater%storage(damID) + flow_vol - flow_res - evap 
     end if

     Trunoff%erout(iunit,nt_nliq) = -flow_res / (theDeltaT)
     !print*, "regulation", damID, StorWater%storage(damID), flow_vol, flow_res, min_flow

     !! evaporation from the reservoir
     !if ( StorWater%storage(damID) > evap ) then
     !  StorWater%storage(damID) = StorWater%storage(damID) - evap
     !else
     !  StorWater%storage(damID) = 0._r8
     !  print*, "WARNING, POTENTIAL EVAP LARGER THEN RESERVOIR STORAGE, iunit ", iunit, evap
     !endif
     ! commented out because evap need to be takren into consideration in the releases
  end subroutine Regulation

!-----------------------------------------------------------------------

  subroutine ExtractionRegulatedFlow(TheDeltaT)
     !! DESCRIPTION: extract water from the reservoir release
     ! !DESCRIPTION: the extraction needs to be distributed accross the dependent unit demand
     !! DESCRIPTION: do not extract more than 10% of the mean monthly flow
     implicit none
     real(r8), intent(in) :: TheDeltaT
     !--- local ---
     integer :: iunit, idam, cnt, ierror, mdam, gdam, begr, endr
     integer :: iter, iflag, iflagm
     integer :: yr, month, day, tod
     logical :: done
     type(mct_aVect) :: aVect_wdG
     real(r8),allocatable :: flow_vol(:)    ! dam storage amount
     real(r8),allocatable :: dam_uptake(:)  ! dam uptake from gridcells locally
     real(r8),allocatable :: dam_uptake_sum(:)  ! dam uptake from gridcells sum over all pes
     real(r8),allocatable :: fracsum(:)     ! sum of dam fraction on gridcells
     real(r8) :: demand, demand_orig, supply
     character(len=*),parameter :: subname='(ExtractionRegulatedFlow)'

     call t_startf('moswrm_ERFlow')

     call get_curr_date(yr, month, day, tod)
     begr = rtmCTL%begr
     endr = rtmCTL%endr

     allocate(dam_uptake(ctlSubwWRM%NDam))
     allocate(dam_uptake_sum(ctlSubwWRM%NDam))
     allocate(fracsum(begr:endr))
     allocate(flow_vol(ctlSubwWRM%LocalNumDam))
     do idam = 1,ctlSubwWRM%LocalNumDam
        iunit = WRMUnit%icell(idam)
        flow_vol(idam) = -Trunoff%erout(iunit,nt_nliq) * theDeltaT
     enddo

!debug     write(iulog,*) subname,' initial flow_vol ',minval(flow_vol),maxval(flow_vol)

     done = .false.
     iter = 0

     !---------------------------
     ! on src side:
     ! avect_wd(1,:) = fraction of demand available
     ! avect_wg(1,:) = current demand from gridcell
     !---------------------------

     do while (.not. done)
        iter = iter + 1

        !--- copy gridcell demand into aVect_wg ---

        call mct_aVect_zero(aVect_wg)
        cnt = 0
        do iunit = begr,endr
           cnt = cnt + 1
           aVect_wg%rAttr(1,cnt) = StorWater%demand(iunit)
        enddo

        !--- sum gridcell demand to dams

        if (barrier_timers) then
           call t_startf('moswrm_ERFlow_avmult_barrier')
           call mpi_barrier(mpicom_rof,ierror)
           call t_stopf('moswrm_ERFlow_avmult_barrier')
        endif

        call mct_aVect_zero(aVect_wd)
        call t_startf('moswrm_ERFlow_avmult')
        call mct_sMat_avMult(aVect_wg, sMatP_g2d, aVect_wd)
        call t_stopf('moswrm_ERFlow_avmult')

        !--- compute new dam fraction from total gridcell demand

!debug        write(iulog,*) subname,' dam demand = ',iam,minval(aVect_wd%rAttr(1,:)),maxval(aVect_wd%rAttr(1,:))
        do idam = 1,ctlSubwWRM%LocalNumDam
           demand = aVect_wd%rAttr(1,idam)
           if (demand > 0._r8) then
!debug              write(iulog,'(2a,2i6,2g20.10)') subname,' volumes ',iam,idam,flow_vol(idam),demand
              aVect_wd%rAttr(1,idam) = flow_vol(idam)/demand
           else
              aVect_wd%rAttr(1,idam) = 0._r8
           endif
        enddo

!debug        write(iulog,'(2a,i8,2g20.10)') subname,' locdam frac =',iam,minval(aVect_wd%rAttr(1,:)),maxval(aVect_wd%rAttr(1,:))

        !--- gather and bcast dam fraction

        if (barrier_timers) then
           call t_startf('moswrm_ERFlow_gather_barrier')
           call mpi_barrier(mpicom_rof,ierror)
           call t_stopf('moswrm_ERFlow_gather_barrier')
        endif

        call t_startf('moswrm_ERFlow_gather')
        call mct_aVect_gather(aVect_wd,aVect_wdG,gsMap_wd,mastertask,mpicom_rof)
        call t_stopf('moswrm_ERFlow_gather')

        if (barrier_timers) then
           call t_startf('moswrm_ERFlow_bcast_barrier')
           call mpi_barrier(mpicom_rof,ierror)
           call t_stopf('moswrm_ERFlow_bcast_barrier')
        endif

        call t_startf('moswrm_ERFlow_bcast')
        call mct_aVect_bcast(aVect_wdG,mastertask,mpicom_rof)
        call t_stopf('moswrm_ERFlow_bcast')

        if (masterproc) then
!debug           write(iulog,'(2a,i8,2g20.10)') subname,' dam frac =',iam,minval(aVect_wdG%rAttr(1,:)),maxval(aVect_wdG%rAttr(1,:))
           !do idam = 1,ctlSubwWRM%NDam
           !   write(iulog,'(2a,2i8,g20.10)') subname,' dam frac =',iam,idam,aVect_wdG%rAttr(1,idam)
           !enddo
        endif

        !--- reduce local demand based on fraction from dam
        !--- if any dam fraction >= 1.0 for a gridcell, then update and iterate
        !--- if any sum of dam fraction >= 1.0 for a gridcell, then update and iterate
        !--- if any sum of dam fraction < 1.0 for a gridcell, then update and iterate

        dam_uptake = 0._r8
        if (maxval(aVect_wdG%rAttr(1,:)) >= 1.0_r8) then
           do iunit = begr, endr
              cnt = 0
              do mdam = 1,WRMUnit%myDamNum(iunit)
                 gdam = WRMUnit%myDam(mdam,iunit)
                 if (aVect_wdG%rAttr(1,gdam) > 1.0_r8) cnt=cnt+1
              enddo
              if (cnt > 0) then
                 demand_orig = StorWater%demand(iunit)
                 do mdam = 1,WRMUnit%myDamNum(iunit)
                    gdam = WRMUnit%myDam(mdam,iunit)
                    if (aVect_wdG%rAttr(1,gdam) > 1.0_r8) then
                       supply = max(demand_orig / float(cnt), StorWater%demand(iunit))
                       dam_uptake(gdam) = dam_uptake(gdam) + supply
                       StorWater%demand(iunit) = StorWater%demand(iunit) - supply
                       StorWater%supply(iunit) = StorWater%supply(iunit) + supply
                    endif
                 enddo
              endif
           enddo
        else
           iflag = 0
           iflagm = 0
           fracsum = 0._r8
           do iunit = begr, endr
              do mdam = 1,WRMUnit%myDamNum(iunit)
                 gdam = WRMUnit%myDam(mdam,iunit)
                 fracsum(iunit) = fracsum(iunit) + aVect_wdG%rAttr(1,gdam)
              enddo
              if (fracsum(iunit) >= 1.0_r8) iflag=1
           enddo
           call shr_mpi_max(iflag,iflagm,mpicom_rof,'wrm iflag',all=.true.)

           do iunit = begr, endr
              if ((iflagm >  0 .and. fracsum(iunit) >= 1.0_r8) .or. &
                  (iflagm == 0 .and. fracsum(iunit) >  0._r8)) then
                 demand_orig = StorWater%demand(iunit)
                 do mdam = 1,WRMUnit%myDamNum(iunit)
                    gdam = WRMUnit%myDam(mdam,iunit)
                    supply = max(demand_orig*aVect_wdG%rAttr(1,gdam)/max(fracsum(iunit),1.0_r8), StorWater%demand(iunit))
                    dam_uptake(gdam) = dam_uptake(gdam) + supply
                    StorWater%demand(iunit) = StorWater%demand(iunit) - supply
                    StorWater%supply(iunit) = StorWater%supply(iunit) + supply
                 enddo
              endif
           enddo
        endif

        call t_startf('moswrm_ERFlow_cleanG')
        call mct_aVect_clean(aVect_wdG)
        call t_stopf('moswrm_ERFlow_cleanG')

        if (barrier_timers) then
           call t_startf('moswrm_ERFlow_sum_barrier')
           call mpi_barrier(mpicom_rof,ierror)
           call t_stopf('moswrm_ERFlow_sum_barrier')
        endif

        call t_startf('moswrm_ERFlow_sum')
        call shr_mpi_sum(dam_uptake,dam_uptake_sum,mpicom_rof,'wrm dam_uptake',all=.true.)
        call t_stopf('moswrm_ERFlow_sum')

        do idam = 1,ctlSubwWRM%LocalNumDam
           gdam = WRMUnit%damID(idam)
           flow_vol(idam) = flow_vol(idam) - dam_uptake_sum(gdam)
        enddo

        if (iter > 2) done = .true.

        if (masterproc) write(iulog,'(2a,i6,g20.10,l4)') subname,' iteration ',iter,sum(dam_uptake_sum),done

     enddo

     do idam = 1,ctlSubwWRM%LocalNumDam
        iunit = WRMUnit%icell(idam)
        Trunoff%erout(iunit,nt_nliq) = - flow_vol(idam) / (theDeltaT)
     enddo
     deallocate(flow_vol)
     deallocate(dam_uptake)
     deallocate(dam_uptake_sum)
     deallocate(fracsum)

     call t_stopf('moswrm_ERFlow')

!------------

#if (1 == 0)

     if (barrier_timers) then
        call t_startf('moswrm_ERFlow_avmult_barrier')
        call mpi_barrier(mpicom_rof,ier)
        call t_stopf('moswrm_ERFlow_avmult_barrier')
     endif

     call t_startf('moswrm_ERFlow_avmult')
     call mct_sMat_avMult(aVect_wg, sMatP_g2d, aVect_wd)
     call t_stopf('moswrm_ERFlow_avmult')

     call t_startf('moswrm_ERFlow_writ1')
     !--- g2d sum ---
     do idam = 1,ctlSubwWRM%LocalNumDam
!        write(iulog,'(2a,2i8,2g20.10)') subname,' idam demand =',idam,WRMUnit%damID(idam),aVect_wd%rAttr(1,idam),StorWater%storage(idam)
     enddo
     call t_stopf('moswrm_ERFlow_writ1')

     !--- copy dam supply into aVect_wd ---
     call t_startf('moswrm_ERFlow_copy2')
     call mct_aVect_zero(aVect_wd)
     do idam = 1,ctlSubwWRM%LocalNumDam
!tcx test        aVect_wd%rAttr(1,cnt) = StorWater%storage(idam)
        aVect_wd%rAttr(1,idam) = WRMUnit%damID(idam)
!        write(iulog,'(2a,i8,g20.10)') subname,' idam setting =',idam,aVect_wd%rAttr(1,idam)
     enddo
     call t_stopf('moswrm_ERFlow_copy2')

     if (barrier_timers) then
        call t_startf('moswrm_ERFlow_gather_barrier')
        call mpi_barrier(mpicom_rof,ier)
        call t_stopf('moswrm_ERFlow_gather_barrier')
     endif

     call t_startf('moswrm_ERFlow_gather')
     call mct_aVect_gather(aVect_wd,aVect_wdG,gsMap_wd,mastertask,mpicom_rof)
     call t_stopf('moswrm_ERFlow_gather')

     if (barrier_timers) then
        call t_startf('moswrm_ERFlow_bcast_barrier')
        call mpi_barrier(mpicom_rof,ier)
        call t_stopf('moswrm_ERFlow_bcast_barrier')
     endif

     call t_startf('moswrm_ERFlow_bcast')
     call mct_aVect_bcast(aVect_wdG,mastertask,mpicom_rof)
     call t_stopf('moswrm_ERFlow_bcast')

     call t_startf('moswrm_ERFlow_writ2')
     !--- dam bcast result ---
     do idam = 1,ctlSubwWRM%NDam
!        write(iulog,'(2a,i8,g20.10)') subname,' idam storage =',idam,aVect_wdG%rAttr(1,idam)
     enddo
     call t_stopf('moswrm_ERFlow_writ2')

     call t_startf('moswrm_ERFlow_cleanG')
     call mct_aVect_clean(aVect_wdG)
     call t_stopf('moswrm_ERFlow_cleanG')

     call t_stopf('moswrm_ERFlow')

!!!!#if (1 == 0)
     demand = 0._r8
     supply = 0._r8 

     damID = WRMUnit%INVicell(iunit)   ! dam index for this gridcell
     if ( damID > ctlSubwWRM%NDam .OR. damID <= 0 ) then
        print*, "Error in Regulation with DamID ",damID
        stop
     end if
     flow_vol = -Trunoff%erout(iunit,nt_nliq) * theDeltaT
     min_flow = 0.1_r8 * WRMUnit%MeanMthFlow(damID, month) * theDeltaT
     !min_flow=10
            
     if ( (flow_vol .le. min_flow) .OR. ( WRMUnit%use_Irrig(damID)  .eq. 0  ) .OR. (flow_vol.lt.0.1_r8) ) then
        !if ( (flow_vol .le. min_flow) ) then
           !print*,"No extraction from regulated flow permitted ", month, damID 
     else

        do idam = 1,WRMUnit%dam_Ndepend(damID)    ! dam dependency counter
           !units where the other reservoirs are located
           ! dam_depend(damID,idam) is the gridcell ID
           ! INVisubw(gridcell_ID) is the local gridcell ID?
           idepend = WRMUnit%INVisubw(WRMUnit%dam_depend(damID,idam))   ! tcx what index are we trying to get here?
!           if ( idepend < 0 .and. idepend <= Tctl%NUnit ) then
           if ( idepend < 0) then
              print*,"ERROR idepend DATABASE extract from reg. flow ",idepend,damID,idam,WRMUnit%dam_depend(damID,idam)
              stop
           end if
           !pro rated demand - no interannual variability
           !demand = demand + ( StorWater%demand(idepend) * WRMUnit%StorCap(damID) / WRMUnit%TotStorCapDepend(idepend))
           !- need adjustement based on storage in each dam
           !if ( StorWater%Storage(damID) > 0._r8 .and. WRMUnit%TotStorCapDepend(idepend) > 0.1_r8 ) then
           !  prorata=1._r8
           !  if ( StorWater%Storage(damID)  / WRMUnit%TotStorCapDepend(idepend) < 1._r8 ) then
           !    prorata = StorWater%Storage(damID)  / WRMUnit%TotStorCapDepend(idepend)
           !  end if
           !  demand = demand +  StorWater%demand(idepend) * prorata
           !end if
!CHANGES HERE FROM HESS2013 PAPER, USE FLOW FOR THER FISTRIBUTION
           prorata=1._r8
           if (WRMUnit%TotInflowDepend(idepend) > 1._r8  .and. WRMUnit%MeanMthFlow(damID, 13) >= 0.001_r8) then
              prorata = WRMUnit%MeanMthFlow(damID, 13)/WRMUnit%TotInflowDepend(idepend)
              if ( prorata > 1._r8) then
                 print*,"prorata >1", prorata,idepend,damID,iunit,WRMUnit%TotInflowDepend(idepend), WRMUnit%MeanMthFlow(damID, 13)
                 prorata = 1._r8
              endif
           end if

           demand = demand +  StorWater%demand(idepend) * prorata
        end do

        if ( (flow_vol - min_flow) >= demand ) then
           supply = demand
        else 
           if ( flow_vol >= min_flow ) then
              supply = flow_vol - min_flow
           end if
        end if
 
        if ( supply .eq. 0._r8  .and. demand > 0._r8) then
           print*, "no regulation, no extraction ",flow_vol, min_flow, supply, demand
        end if
        !if ( supply .lt. 0 ) then
        !  print*,"Error computing extraction from reservoir"
        !  print*, supply, demand, min_flow, flow_vol
        !  stop
        !end if
                  
! SECOND LOOP
        do idam = 1,WRMUnit%dam_Ndepend(damID)
           idepend = WRMUnit%INVisubw(WRMUnit%dam_depend(damID,idam))
           ! pro rated supply
           if ( idepend > 0 .and. demand > 0._r8 .and. WRMUnit%TotStorCapDepend(idepend) > 0.1_r8 .and. WRMUnit%MeanMthFlow(damID, 13) >= 0.001_r8)  then 
              !StorWater%demand(idepend)= StorWater%demand(idepend) - (supply/demand * StorWater%demand(idepend) / WRMUnit%subw_Ndepend(idepend))
              !StorWater%supply(idepend)= StorWater%supply(idepend) + (supply/demand * StorWater%demand(idepend) / WRMUnit%subw_Ndepend(idepend))
              !StorWater%demand(idepend)= StorWater%demand(idepend) - (supply/demand * StorWater%demand(idepend) * WRMUnit%StorCap(damID) / WRMUnit%TotStorCapDepend(idepend))
              !StorWater%supply(idepend)= StorWater%supply(idepend) + (supply/demand * StorWater%demand(idepend) * WRMUnit%StorCap(damID) / WRMUnit%TotStorCapDepend(idepend))
              !CHANGES start
              !prorata=1._r8
              !if ( StorWater%Storage(damID)  / WRMUnit%TotStorCapDepend(idepend) < 1._r8 ) then
              !  prorata = StorWater%Storage(damID)  / WRMUnit%TotStorCapDepend(idepend)
              !end if
              !if ( prorata < 0._r8 ) prorata = 0._r8
!CHANGES HERE FROM HESS2013 PAPER, USE FLOW FOR THER FISTRIBUTION
              prorata = 1._r8
              if (WRMUnit%TotInflowDepend(idepend) > 1._r8  ) then
                 prorata = WRMUnit%MeanMthFlow(damID, 13)/WRMUnit%TotInflowDepend(idepend)
                 if (prorata > 1._r8) then
                    print*, "Error prorata ",prorata, idepend, WRMUnit%TotInflowDepend(idepend), WRMUnit%MeanMthFlow(damID, 13)
                    prorata = 1._r8
                 endif
              end if

              StorWater%supply(idepend)= StorWater%supply(idepend) + (supply/demand * StorWater%demand(idepend) * prorata)
              StorWater%demand(idepend)= StorWater%demand(idepend) - (supply/demand * StorWater%demand(idepend) * prorata)
              if ( StorWater%supply(idepend) .lt. 0._r8 ) then
                 print*, "error supply"
                 stop
              else if ( StorWater%demand(idepend) .lt. 0._r8 ) then
                 print*,"Error demand in first loop", StorWater%demand(idepend), supply/demand
                 stop
              end if  
           end if
        end do

        flow_vol = flow_vol - supply
        Trunoff%erout(iunit,nt_nliq) = - flow_vol / (theDeltaT)

        supply = 0._r8
        demand = 0._r8

        do idam = 1,WRMUnit%dam_Ndepend(damID)
           idepend = WRMUnit%INVisubw(WRMUnit%dam_depend(damID,idam))
           if ( StorWater%demand(idepend) > 0._r8) then
              demand = demand + StorWater%demand(idepend) ! total demand with no prorata
           endif
        end do

        if ( demand > 0._r8) then
           if ( (flow_vol - min_flow) >= demand ) then
              supply = demand
           else
              if ( flow_vol > min_flow ) then
                 supply = flow_vol - min_flow
              end if
           end if

           do idam = 1,WRMUnit%dam_Ndepend(damID)
              idepend = WRMUnit%INVisubw(WRMUnit%dam_depend(damID,idam))
              ! pro rated supply
              if ( idepend > 0 .and. StorWater%demand(idepend) > 0._r8 ) then
                 StorWater%supply(idepend)= StorWater%supply(idepend) + (supply/demand * StorWater%demand(idepend))
                 !StorWater%demand(idepend)= StorWater%demand(idepend) - (supply/demand * StorWater%demand(idepend))
                 StorWater%demand(idepend)= StorWater%demand(idepend) * ( 1._r8 - supply/demand )
                 if ( StorWater%demand(idepend) .lt. 0._r8 ) then
                    print*,"Error demand", StorWater%demand(idepend), supply/demand, supply, demand
                    stop
                 end if
              end if
           end do
        endif

! END SECOND LOOP
        flow_vol = flow_vol - supply
        Trunoff%erout(iunit,nt_nliq) = - flow_vol / (theDeltaT)
     end if
#endif

  end subroutine ExtractionRegulatedFlow

!-----------------------------------------------------------------------
end MODULE WRM_modules
