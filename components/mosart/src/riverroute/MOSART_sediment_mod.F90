!
MODULE MOSART_sediment_mod
! Description: core code of MOSART-sediment. Can be incoporated within any land model via a interface module
! 
! Developed by Hongyi Li, 03/2015.
! REVISION HISTORY:
!-----------------------------------------------------------------------

! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
    use shr_const_mod , only : denh2o => SHR_CONST_RHOFW, denice => SHR_CONST_RHOICE, grav => SHR_CONST_G, SHR_CONST_REARTH, SHR_CONST_PI
    use RtmVar        , only : iulog, inst_suffix, smat_option
    use shr_sys_mod   , only : shr_sys_abort
    !use clm_varcon , only : denh2o, grav !!density of liquid water [kg/m3], gravity constant [m/s2]
    use RunoffMod, only : Tctl, TUnit, TRunoff, TPara
    use RunoffMod, only : rtmCTL
    use rof_cpl_indices, only : nt_rtm, rtm_tracers, nt_nliq, nt_nice, nt_nmud, nt_nsan, KW, DW
    use MOSART_BGC_type, only : TSedi, TSedi_para
    implicit none
    real(r8), parameter :: TINYVALUE_s = 1.0e-12_r8  ! double precision variable has a significance of about 16 decimal digits
    real(r8), parameter :: densedi = 2650._r8      ! density of sediment [Kg/m3]
    !real(r8), parameter :: diam_sand = 0.00035_r8  ! grain size of sand, 0.35mm
    real(r8), parameter :: KVIS = 1.0036e-6_r8     ! the kinematic viscosity (m2/s)at 20°Ê*/
    real(r8), parameter :: DVIS = 1.002e-3_r8      ! the dynamatic viscosity (Pa°§s)at 20°Ê*/

! !PUBLIC MEMBER FUNCTIONS:
    contains

    subroutine hillslopeSediment(iunit, theDeltaT)
    ! !DESCRIPTION: Hillslope erosion/routing of sediment
        implicit none

        integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT

        real(r8) :: Tau_h
        TRunoff%ehout(iunit, nt_nmud) = -TRunoff%qsur(iunit,nt_nmud) !0._r8
        TRunoff%qsub(iunit, nt_nmud) = 0._r8
        !TRunoff%etin(iunit,nt_nmud) = 0._r8
        !TRunoff%etin(iunit,nt_nmud) = TRunoff%qsur(iunit,nt_nmud) * TUnit%area(iunit)

        !TRunoff%etin(iunit,nt_nmud) = 1.e-6_r8 * TUnit%area(iunit) * TUnit%frac(iunit) ! Ks/s
        !TRunoff%etin(iunit,nt_nsan) = 1.e-6_r8 * TUnit%area(iunit) * TUnit%frac(iunit) ! Ks/s
    end subroutine hillslopeSediment

    subroutine subnetworkSediment(iunit, theDeltaT)
    ! !DESCRIPTION: subnetwork channel processes for mud sediment. Assume no sand sediment in tributary channels
        implicit none    
        integer, intent(in) :: iunit !, nt      
        real(r8), intent(in) :: theDeltaT

        !real(r8) :: ers       ! sand-sediment erosion
        !real(r8) :: ersal     ! sand-sediment erosion @ active layer
        !real(r8) :: ersb      ! sand-sediment erosion @ bank
        !real(r8) :: ses       ! sand-sediment settled to bed zone
        !real(r8) :: seout     ! sand-sediment discharge out of local channel
        real(r8) :: erm       ! mud-sediment erosion
        real(r8) :: ermal     ! mud-sediment erosion @ active layer
        real(r8) :: ermb      ! mud-sediment erosion @ bank
        real(r8) :: sem       ! mud-sediment settled to bed zone
        real(r8) :: meout     ! mud-sediment discharge out of local channel
        real(r8) :: erb       ! total-sediment erosion @ bank
        real(r8) :: temp1, temp2, w_temp


        real(r8) :: Achannel  ! effective channel area

            ! sand-sediment
            !ses = CRSES(iunit,TRunoff%yt(iunit,nt_nliq),TUnit%tslp(iunit),TRunoff%conc_t(iunit,nt_nsan))
            !ers = CRERS(TRunoff%yt(iunit,nt_nliq),TUnit%tslp(iunit))
            !seout = -TRunoff%etout(iunit,nt_nliq)*TRunoff%conc_t(iunit,nt_nsan)
            ! mud-sediment
            sem = 0._r8 !CRSEM(TRunoff%yt(iunit,nt_nliq),TUnit%tslp(iunit),TRunoff%conc_t(iunit,nt_nmud))
            erm = 0._r8 !CRERM(TRunoff%yt(iunit,nt_nliq),TUnit%tslp(iunit),TRunoff%conc_t(iunit,nt_nmud))
            meout = -TRunoff%etout(iunit,nt_nliq)*TRunoff%conc_t(iunit,nt_nmud)

            if(TRunoff%wt(iunit,nt_nmud) <= TINYVALUE_s) then
                meout = 0._r8
                sem = 0._r8
            else if ((meout+sem) * theDeltaT > TRunoff%wt(iunit,nt_nmud)+TINYVALUE_s) then
                w_temp = TRunoff%wt(iunit,nt_nmud) * 0.95_r8
                temp1 = sem/(sem+meout)
                temp2 = meout/(sem+meout)
                sem = temp1*w_temp/theDeltaT
                meout = temp2*w_temp/theDeltaT
            end if

            erm = 0._r8
            ermal = 0._r8
            ermb = 0._r8
            sem = 0._r8
            TSedi%ermal_t(iunit)  =   ermal
            TSedi%sem_t(iunit)    =   sem
            TSedi%ermb_t(iunit)   =   ermb            

            TRunoff%etout(iunit,nt_nsan) = 0._r8
            TRunoff%etout(iunit,nt_nmud) = -meout

            TSedi%Ssal_t(iunit) = 0._r8 ! TSedi%Ssal_t(iunit) + (TSedi%ses_t(iunit) - TSedi%ersal_t(iunit)) * theDeltaT
            !TSedi%Smal_t(iunit) = 0._r8 ! TSedi%Smal_t(iunit) + (TSedi%sem_t(iunit) - TSedi%ermal_t(iunit)) * theDeltaT

            TRunoff%dwt(iunit,nt_nsan) = 0._r8 !TRunoff%etin(iunit,nt_nsan) + TRunoff%etout(iunit,nt_nsan) + TSedi%ersal_t(iunit) - TSedi%ses_t(iunit) + TSedi%ersb_t(iunit)
            TRunoff%dwt(iunit,nt_nmud) = TRunoff%etin(iunit,nt_nmud) + TRunoff%etout(iunit,nt_nmud) + TSedi%ermal_t(iunit) - TSedi%sem_t(iunit) + TSedi%ermb_t(iunit)

            !TRunoff%wt(iunit,nt_nsan) = TRunoff%wt(iunit,nt_nsan) + TRunoff%dwt(iunit,nt_nsan)*theDeltaT
            !TRunoff%wt(iunit,nt_nmud) = TRunoff%wt(iunit,nt_nmud) + TRunoff%dwt(iunit,nt_nmud)*theDeltaT

            !! for budget calculation
            TRunoff%wt_al(iunit,nt_nsan) = 0._r8 ! TSedi%Ssal_t(iunit)
            !TRunoff%wt_al(iunit,nt_nmud) = TSedi%Smal_t(iunit)
            !TRunoff%etexchange(iunit,nt_nsan) = TSedi%ersb_t(iunit)
            !TRunoff%etexchange(iunit,nt_nmud) = TSedi%ermb_t(iunit)

            if(TRunoff%wt(iunit,nt_nmud) < -1.e-10) then
               write(iulog,*) 'Negative mud storage in t-zone! ', iunit, TRunoff%wt(iunit, nt_nmud) , TRunoff%etin(iunit,nt_nmud), TRunoff%etout(iunit,nt_nmud), TSedi%ermal_t(iunit), TSedi%sem_t(iunit), TSedi%ermb_t(iunit)
               call shr_sys_abort('mosart: negative mud t-zone storage')
            end if
    end subroutine subnetworkSediment

    subroutine mainchannelSediment(iunit, theDeltaT)
    ! !DESCRIPTION: subnetwork channel routing
        implicit none    
        integer, intent(in) :: iunit !,nt      
        real(r8), intent(in) :: theDeltaT

        real(r8) :: ers       ! sand-sediment erosion
        real(r8) :: ersal     ! sand-sediment erosion @ active layer
        real(r8) :: ersb      ! sand-sediment erosion @ bank
        real(r8) :: ses       ! sand-sediment settled to bed zone
        real(r8) :: seout     ! sand-sediment discharge out of local channel
        real(r8) :: erm       ! mud-sediment erosion
        real(r8) :: ermal     ! mud-sediment erosion @ active layer
        real(r8) :: ermb      ! mud-sediment erosion @ bank
        real(r8) :: sem       ! mud-sediment settled to bed zone
        real(r8) :: meout     ! mud-sediment discharge out of local channel
        real(r8) :: erb       ! total-sediment erosion @ bank
        real(r8) :: ers_net   ! net erosion rate
        real(r8) :: temp1, temp2, w_temp

        real(r8) :: Achannel  ! effective channel area

        integer :: k, nt, myflag
        real(r8) :: s_conc_equi! sand-sediment concentration under the equilibirum state (erosion rate = deposition rate)
        real(r8) :: s_wr_equi! sand-sediment storage under the equilibirum state (erosion rate = deposition rate)

        TRunoff%erin(iunit,nt_nsan) = -TRunoff%eroutUp(iunit,nt_nsan)
        TRunoff%erin(iunit,nt_nmud) = -TRunoff%eroutUp(iunit,nt_nmud)

        ! bed-material load (sand or coarser)
        if(TRunoff%vr(iunit,nt_nliq)>TINYVALUE_s) then
            seout = CRQs_EH(TRunoff%yr(iunit,nt_nliq), abs(TRunoff%vr(iunit,nt_nliq)), TUnit%rwidth(iunit), TUnit%rslp(iunit), TUnit%nr(iunit), TSedi_para%d50(iunit))
            !seout = CRQs_Wu(TRunoff%yr(iunit,nt_nliq), TRunoff%vr(iunit,nt_nliq), TUnit%rwidth(iunit), TUnit%rslp(iunit), TSedi_para%d50(iunit))
            !seout = CRQs_Parker(TRunoff%yr(iunit,nt_nliq), TRunoff%vr(iunit,nt_nliq), TUnit%rwidth(iunit), TUnit%rslp(iunit), TSedi_para%d50(iunit))
         else
            seout = 0._r8
        end if        

        ! wash load (mud or finer)
        if(TRunoff%wr(iunit,nt_nliq)<TINYVALUE_s .and. abs(TRunoff%erout(iunit,nt_nliq))>TINYVALUE_s) then
            meout = TRunoff%erin(iunit,nt_nmud) + TRunoff%erlateral(iunit,nt_nmud)
        else
            meout = -TRunoff%erout(iunit,nt_nliq)*TRunoff%conc_r(iunit,nt_nmud)     
            if (meout * theDeltaT > TRunoff%wr(iunit,nt_nmud)+TINYVALUE_s) then
                meout = TRunoff%wr(iunit,nt_nmud)*0.95_r8/theDeltaT
            end if        
        end if

        TRunoff%erout(iunit,nt_nsan) = -seout
        TRunoff%erout(iunit,nt_nmud) = -meout
        ! in case diffusion wave method is used, account for the mass balance caused by backwater effects

        !== Please don't remove the block
        !==TODO: in current MOSART, the upstream/downstream relationship and mass balance are too difficult to track
        !        hence the block below is causing some mass balance issue

        !if(Tctl%RoutingMethod == 4) then
        !do nt=nt_nmud,nt_nsan
        !
        !    if(TRunoff%rslp_energy(iunit) >= TINYVALUE_s) then ! flow is from current channel to downstream
        !      if(TRunoff%erin(iunit,nt)*theDeltaT + TRunoff%wr(iunit,nt) <= TINYVALUE_s) then! too much negative inflow from upstream, 
        !         TRunoff%erout(iunit,nt) = 0._r8
        !      elseif(TRunoff%erout(iunit,nt) <= -TINYVALUE_s .and. TRunoff%wr(iunit,nt) + &
        !         (TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT < TINYVALUE_s) then
        !         TRunoff%erout(iunit,nt) = -(TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%wr(iunit,nt)*0.95_r8 / theDeltaT)
        !      end if
        !    elseif(TRunoff%rslp_energy(iunit) <= -TINYVALUE_s) then ! flow is from downstream to current channel
        !       if(rtmCTL%nUp_dstrm(iunit) > 1) then
        !           if(TRunoff%erin_dstrm(iunit,nt)*theDeltaT + TRunoff%wr_dstrm(iunit,nt)/rtmCTL%nUp_dstrm(iunit) <= TINYVALUE_s) then! much negative inflow from upstream,
        !               TRunoff%erout(iunit,nt) = 0._r8
        !          elseif(TRunoff%erout(iunit,nt) >= TINYVALUE_s .and. TRunoff%wr_dstrm(iunit,nt)/rtmCTL%nUp_dstrm(iunit) &
        !              - TRunoff%erout(iunit,nt) * theDeltaT < TINYVALUE_s) then
        !              TRunoff%erout(iunit,nt) = TRunoff%wr_dstrm(iunit,nt)*0.95_r8 / theDeltaT / rtmCTL%nUp_dstrm(iunit)
        !           end if    
        !       else
        !           if(TRunoff%erin_dstrm(iunit,nt)*theDeltaT + TRunoff%wr_dstrm(iunit,nt) <= TINYVALUE_s) then! much negative inflow from upstream,
        !               TRunoff%erout(iunit,nt) = 0._r8
        !           elseif(TRunoff%erout(iunit,nt) >= TINYVALUE_s .and. TRunoff%wr_dstrm(iunit,nt) &
        !             - TRunoff%erout(iunit,nt) * theDeltaT < TINYVALUE_s) then
        !              TRunoff%erout(iunit,nt) = TRunoff%wr_dstrm(iunit,nt)*0.95_r8 / theDeltaT
        !           end if
        !       end if                 
        !       
        !    else  ! no flow between current channel and downstream
        !      TRunoff%erout(iunit,nt) = 0._r8
        !    end if
        !
        !end do
        !end if
        !== Please don't remove the block

        !==TODO: in current MOSART, the upstream/downstream relationship and mass balance are too difficult to track, hence assuming there is no sediment flux from downstream to upstream channel when backwater occurs
        if(Tctl%RoutingMethod == DW) then
        do nt=nt_nmud,nt_nsan

            if(abs(TRunoff%rslp_energy(iunit)) < TINYVALUE_s) then ! no flow between current channel and downstream
              TRunoff%erout(iunit,nt) = 0._r8
            elseif(TRunoff%rslp_energy(iunit) <= -TINYVALUE_s) then ! flow is from downstream to current channel
              TRunoff%erout(iunit,nt) = 0._r8                
            end if
        end do
        end if

        seout = - TRunoff%erout(iunit,nt_nsan) 
        meout = - TRunoff%erout(iunit,nt_nmud) 

        if(abs(TRunoff%erout(iunit,nt_nliq)) > TINYVALUE_s) then
            s_conc_equi = seout/abs(TRunoff%erout(iunit,nt_nliq))        
            s_wr_equi = s_conc_equi * TRunoff%wr(iunit,nt_nliq)
        else 
            s_wr_equi = 0._r8
        end if

        ers = 0._r8
        ses = 0._r8
        if(s_wr_equi > TRunoff%wr(iunit,nt_nsan) + TINYVALUE_s) then ! too less sand in the water, net erosion occurs
            ers = (s_wr_equi - TRunoff%wr(iunit,nt_nsan))/theDeltaT
            ses = 0._r8
        else
            if(s_wr_equi < TRunoff%wr(iunit,nt_nsan) - TINYVALUE_s) then ! too much sand in the water, net deposition occurs
                ers = 0._r8
                ses = (TRunoff%wr(iunit,nt_nsan) - s_wr_equi)/theDeltaT
            end if
        end if

        ! Don't change seout here, only adjust ses and ers for mass balance
        if ((seout+ses) * theDeltaT > (TRunoff%wr(iunit,nt_nsan) + ers* theDeltaT + TINYVALUE_s)) then
            if (seout * theDeltaT < (TRunoff%wr(iunit,nt_nsan) + ers* theDeltaT - TINYVALUE_s)) then
                ses = TRunoff%wr(iunit,nt_nsan)/theDeltaT + ers - seout
                ers = 0._r8
            end if
            if (seout * theDeltaT > (TRunoff%wr(iunit,nt_nsan) + ers* theDeltaT + TINYVALUE_s)) then
                ses = 0._r8
                ers = seout - TRunoff%wr(iunit,nt_nsan)/theDeltaT
            end if
        end if

        if(TRunoff%wr(iunit,nt_nsan) <= TINYVALUE_s) then
            ses = 0._r8             
        end if

        TSedi%Ssal_r(iunit) = TRunoff%wr_al(iunit, nt_nsan)
        if(TSedi%Ssal_r(iunit) >= ers*theDeltaT + TINYVALUE_s) then ! first erosion from active layer, it no sufficient storage, from channel bed
            ersal = ers
        else
            ersal = TSedi%Ssal_r(iunit)*0.95_r8/theDeltaT          
        end if
        ersb = ers - ersal

         ! mud-sediment
        sem = 0._r8 !CRSEM(TRunoff%yr(iunit,nt_nliq),TUnit%rslp(iunit),TRunoff%conc_r(iunit,nt_nmud))
        erm = 0._r8 !CRERM(TRunoff%yr(iunit,nt_nliq),TUnit%rslp(iunit),TRunoff%conc_r(iunit,nt_nmud))
        ermb = 0._r8
        ermal = 0._r8

        TSedi%ers_r(iunit) = ers
        TSedi%ersal_r(iunit) = ersal
        TSedi%ersb_r(iunit) = ersb
        TSedi%ses_r(iunit) = ses
        TSedi%erm_r(iunit) = erm
        TSedi%ermal_r(iunit) = ermal
        TSedi%ermb_r(iunit) = ermb
        TSedi%sem_r(iunit) = sem

        TRunoff%dwr_al(iunit,nt_nsan) = TSedi%ses_r(iunit) - TSedi%ersal_r(iunit)
        TRunoff%dwr_al(iunit,nt_nmud) = 0._r8 !TSedi%sem_r(iunit) - TSedi%ermal_r(iunit)

        !TSedi%Ssal_r(iunit) = TSedi%Ssal_r(iunit) + TRunoff%dwr_al(iunit,nt_nsan) * theDeltaT
        !TSedi%Smal_r(iunit) = TSedi%Smal_r(iunit) + (TSedi%sem_r(iunit) - TSedi%ermal_r(iunit)) * theDeltaT

        TRunoff%dwr(iunit,nt_nsan) = TRunoff%erlateral(iunit,nt_nsan) + TRunoff%erin(iunit,nt_nsan) + TRunoff%erout(iunit,nt_nsan) + TSedi%ersal_r(iunit) - TSedi%ses_r(iunit) + TSedi%ersb_r(iunit)
        TRunoff%dwr(iunit,nt_nmud) = TRunoff%erlateral(iunit,nt_nmud) + TRunoff%erin(iunit,nt_nmud) + TRunoff%erout(iunit,nt_nmud) + TSedi%ermal_r(iunit) - TSedi%sem_r(iunit) + TSedi%ermb_r(iunit)

        if(TRunoff%wr(iunit,nt_nmud).gt.TINYVALUE_s .and. (TRunoff%wr(iunit,nt_nmud) + TRunoff%dwr(iunit,nt_nmud)*theDeltaT)/TRunoff%wr(iunit,nt_nmud) < -1.e-8 .and. (TRunoff%wr(iunit,nt_nmud) + TRunoff%dwr(iunit,nt_nmud)*theDeltaT) < -1.e-8) then
           write(iulog,*) 'Negative mud storage in r-zone! ', iunit, TRunoff%wr(iunit, nt_nmud), TRunoff%erlateral(iunit,nt_nmud), TRunoff%erin(iunit,nt_nmud), TRunoff%erout(iunit,nt_nmud), TSedi%ermal_r(iunit), - TSedi%sem_r(iunit), TSedi%ermb_r(iunit) 
            !write(unit=5112,fmt="((i10), 7(e12.3))") iunit, -TRunoff%erout(iunit,nt_nliq), TRunoff%yr(iunit,nt_nliq), abs(TRunoff%vr(iunit,nt_nliq)), TUnit%rwidth(iunit), TUnit%rslp(iunit), TUnit%nr(iunit), TSedi_para%d50(iunit)
            !write(unit=5113,fmt="((i10), 7(e12.3))") iunit, TRunoff%wr(iunit,nt_nmud), meout, sem, erm, ermal, ermb
           call shr_sys_abort('mosart: negative mud r-zone storage')
        end if
        if(TRunoff%wr(iunit,nt_nsan).gt.TINYVALUE_s .and. (TRunoff%wr(iunit,nt_nsan) + TRunoff%dwr(iunit,nt_nsan)*theDeltaT)/TRunoff%wr(iunit,nt_nsan) < -1.e-8 .and. (TRunoff%wr(iunit,nt_nsan) + TRunoff%dwr(iunit,nt_nsan)*theDeltaT) < -1.e-8) then
           write(iulog,*) 'Negative sand storage in r-zone! ', iunit, TRunoff%wr(iunit, nt_nsan), TRunoff%erlateral(iunit,nt_nsan), TRunoff%erin(iunit,nt_nsan), TRunoff%erout(iunit,nt_nsan), TSedi%ersal_r(iunit), - TSedi%ses_r(iunit), TSedi%ersb_r(iunit) 
            !write(unit=1110,fmt="(i10, 6(e12.3))") iunit,TRunoff%erlateral(iunit,nt_nsan), TRunoff%erin(iunit,nt_nsan), TRunoff%erout(iunit,nt_nsan), TSedi%ersal_r(iunit), - TSedi%ses_r(iunit), TSedi%ersb_r(iunit)
            !write(unit=6112,fmt="((i10), 7(e12.3))") iunit, -TRunoff%erout(iunit,nt_nliq), TRunoff%yr(iunit,nt_nliq), abs(TRunoff%vr(iunit,nt_nliq)), TUnit%rwidth(iunit), TUnit%rslp(iunit), TUnit%nr(iunit), TSedi_para%d50(iunit)
            !write(unit=6113,fmt="((i10), 7(e12.3))") iunit, TRunoff%wr(iunit,nt_nsan), s_wr_equi, seout, ses, ers, ersal, ersb
            !write(unit=6114,fmt="((i10), 5(e12.3))") iunit, TRunoff%wr(iunit,nt_nsan) + TRunoff%dwr(iunit,nt_nsan)*theDeltaT, TRunoff%wr(iunit,nt_nsan), s_wr_equi, (TRunoff%wr(iunit,nt_nsan) - s_wr_equi)/theDeltaT, ses
           call shr_sys_abort('mosart: negative sand r-zone storage')
        end if

        if(TSedi%Ssal_r(iunit) < -1.e-10) then
           write(iulog,*) 'Negative storage of sand active layer  in r-zone ! ', iunit, TSedi%Ssal_r(iunit), s_conc_equi, TRunoff%wr(iunit,nt_nliq), TRunoff%wr(iunit,nt_nsan)
           call shr_sys_abort('mosart: negative sand r-zone active layer storage')
        end if

            !if((TRunoff%yr(iunit,nt_nliq).gt.TINYVALUE_s) .and. (s_wr_equi.gt.TINYVALUE_s)) then
            !    write(unit=1110,fmt="(i10, 6(e12.3))") iunit, ers, ses, seout, TRunoff%wr(iunit,nt_nsan), TRunoff%vr(iunit,nt_nliq), CRQs(TRunoff%yr(iunit,nt_nliq), TRunoff%vr(iunit,nt_nliq), TUnit%rwidth(iunit), TUnit%rslp(iunit), TUnit%nr(iunit))
            !end if

            !if((TRunoff%yr(iunit,nt_nliq).gt.TINYVALUE_s) .and. (TRunoff%wr(iunit,nt_nsan).gt.TINYVALUE_s)) then
            !    write(unit=1111,fmt="(i10, 6(e12.3))") iunit, ers, ses, seout, TRunoff%wr(iunit,nt_nsan), TRunoff%vr(iunit,nt_nliq), CRQs(TRunoff%yr(iunit,nt_nliq), TRunoff%vr(iunit,nt_nliq), TUnit%rwidth(iunit), TUnit%rslp(iunit), TUnit%nr(iunit))
            !end if

    end subroutine mainchannelSediment

    function CRQs_Wu(h_, U_, rwidth_,slope_, D50_) result(Qs_)
      ! !DESCRIPTION: calculate equilibrium sand-sediment transport rate using Wu (2000) equation
      implicit none
      real(r8), intent(in)  :: h_    !channel water depth [m]
      real(r8), intent(in)  :: U_    !channel velocity [m/s]
      real(r8), intent(in)  :: rwidth_!bankfull width [m]
      real(r8), intent(in)  :: slope_  ! channel slope [-]
      real(r8), intent(in)  :: D50_  ! median sediment particle size [m]

      real(r8)             :: Qs_       ! equilibrium sediment transporate rate, [kg/s]

      real(r8) :: R_     ! submerged specific gravity of sediment [-]
      real(r8) :: vsf_    ! sediment settling velocity [m/s]
      real(r8) :: tah_b  !the bed shear stress
      real(r8) :: theta ! Shields parameter [-]
      real(r8) :: theta_c ! critical shear stress [-]
      real(r8) :: phi_   ! dimensionless shear stress [-]
      real(r8) :: phi_s ! dimensionless sediment load
      real(r8) :: q_s    ! the sediment volumetric discharge per unit width, [m2/s] 

      R_ = (densedi-denh2o)/denh2o
      theta = h_*slope_/((R_ - 1._r8)*D50_)
      theta_c = 0.0386_r8
      phi_ = theta/theta_c

      if(phi_ <= 1._r8) then
          Qs_ = 0._r8
          return
      end if      
      vsf_ = CRVSF(D50_)
      phi_s = 0.0000262_r8 * ((phi_ - 1._r8)*U_/vsf_)**1.74_r8

      q_s = phi_s * sqrt(grav*R_*D50_)*D50_
      Qs_ = densedi * rwidth_ * q_s

      if(abs(Qs_) < TINYVALUE_s) then
          Qs_ = 0._r8
      endif

      return
    end function CRQs_Wu    


    function CRQs_Parker(h_, U_, rwidth_,slope_, D50_) result(Qs_)
      ! !DESCRIPTION: calculate equilibrium sand-sediment transport rate using Parker (1990) equation
      implicit none
      real(r8), intent(in)  :: h_    !channel water depth [m]
      real(r8), intent(in)  :: U_    !channel velocity [m/s]
      real(r8), intent(in)  :: rwidth_!bankfull width [m]
      real(r8), intent(in)  :: slope_  ! channel slope [-]
      real(r8), intent(in)  :: D50_  ! median sediment particle size [m]
      real(r8)             :: Qs_       ! equilibrium sediment transporate rate, [kg/s]

      real(r8) :: R_     ! submerged specific gravity of sediment [-]
      real(r8) :: tah_b  !the bed shear stress
      real(r8) :: theta ! Shields parameter [-]
      real(r8) :: theta_c ! critical shear stress [-]
      real(r8) :: phi_ ! critical shear stress [-]
      real(r8) :: q_s    ! the sediment volumetric discharge per unit width, [m2/s] 
      real(r8) :: susp_ratio   ! ratio of suspended sediment load over the total sediment load [-]

      R_ = (densedi-denh2o)/denh2o
      theta = h_*slope_/((R_ - 1._r8)*D50_)
      theta_c = 0.0386_r8
      phi_ = theta/theta_c

      q_s = CRf_G(phi_) *(grav*h_*slope_)**1.5/(grav*(R_-1._r8))
      Qs_ = densedi * rwidth_ * q_s

      if(abs(Qs_) < TINYVALUE_s) then
          Qs_ = 0._r8
      endif
      susp_ratio = CRsuspended_ratio(D50_, h_, slope_)
      Qs_ = Qs_ * susp_ratio

      return
    end function CRQs_Parker    

    function CRf_G(phi_) result(G_)
      ! !DESCRIPTION: part of the Parker (1990) equation for total sediment load
      implicit none
      real(r8), intent(in) :: phi_  ! [-]
      real(r8)             :: G_       ! 

      if(phi_ > 1.59_r8) then
          G_ = 11.933_r8*(1._r8 - 0.853_r8/phi_)**4.5_r8
      else
          if(phi_ < 1.0) then
              G_ = 0.00218_r8*phi_**14.2_r8 
          else
              G_ = 0.00218_r8*exp(14.2_r8*(phi_-1._r8) - 9.28_r8*(phi_ - 1._r8)**2._r8)
          end if
      end if

      return
    end function CRf_G

    function CRQs_EH(h_, U_, rwidth_,slope_, roughness_, D50_) result(Qs_)
      ! !DESCRIPTION: calculate equilibrium sand-sediment transport rate using classic Engelund-Hansen equation
      implicit none
      real(r8), intent(in)  :: h_    !channel water depth [m]
      real(r8), intent(in)  :: U_    !channel velocity [m/s]
      real(r8), intent(in)  :: rwidth_!bankfull width [m]
      real(r8), intent(in)  :: slope_  ! channel slope [-]
      real(r8), intent(in)  :: roughness_  ! Manning's roughness [-]
      real(r8), intent(in)  :: D50_  ! median sediment particle size [m]
      real(r8)              :: Qs_       ! equilibrium sediment transporate rate, [kg/s]

      real(r8) :: Cf     ! resistence coefficient
      real(r8) :: R_     ! submerged specific gravity of sediment [-]
      real(r8) :: tah_star !the Shields number (dimensionless shear stress)
      real(r8) :: q_star ! Einstein number (dimensionless sediment flux)    
      real(r8) :: q_s    ! equilibrium the sediment volumetric discharge per unit width, [m2/s] 
      real(r8) :: susp_ratio   ! ratio of suspended sediment load over the total sediment load [-]

      R_ = (densedi-denh2o)/denh2o
      Cf = grav*(roughness_**2)/(h_**(1._r8/3._r8))
      tah_star = Cf * U_* U_ /(R_*grav*D50_)
      q_star  = 0.05_r8 * (tah_star**2.5_r8)/Cf      
      q_s = q_star * sqrt(R_*grav*D50_)*D50_  
      Qs_ = densedi * rwidth_ * q_s

      if(abs(Qs_) < TINYVALUE_s) then
          Qs_ = 0._r8
      endif
      susp_ratio = CRsuspended_ratio(D50_, h_, slope_)
      Qs_ = Qs_ * susp_ratio
            !if(h_>TINYVALUE_s .and. U_>TINYVALUE_s) then
            !    write(unit=1112,fmt="(5(e20.11))") Cf, tah_star, q_star, q_s, Qs_
            !end if
      return
    end function CRQs_EH 


    function CRsuspended_ratio(D_, h_, slope_) result(ratio_)
      ! !DESCRIPTION: calculate the ratio of suspended portion to the total sediment load
      implicit none
      real(r8), intent(in) :: h_,slope_  ! channel water depth [m], channel slope [-]
      real(r8), intent(in) :: D_   !  particle size [m]
      real(r8)             :: ratio_       ! ratio of suspended portion to the total sediment load [-]

      real(r8) :: u_t     ! friction shear velocity [m/s]
      real(r8) :: vsf_    ! sediment settling velocity [m/s]
      real(r8) :: k_      ! von Kármán constant [-]
      real(r8) :: Z_      ! [-]

      k_ = 0.41_r8
      vsf_ = CRVSF(D_)
      u_t = sqrt(grav*h_*slope_)
      if(u_t < TINYVALUE_s) then
          ratio_ = 1._r8
      else
          Z_ = vsf_/(k_*u_t)
          ratio_ = 2.5_r8 * exp(-Z_)
      end if
      if(abs(ratio_) < TINYVALUE_s) then
          ratio_ = 0._r8
      endif
      if(abs(ratio_) > 1._r8) then
          ratio_ = 1._r8
      endif

      return
    end function CRsuspended_ratio    


    function CRVSF(D_) result(vsf_)
      ! !DESCRIPTION: calculate particle settling velocity [m/s] using the equation from Cheng N. (1997). Simplified settling velocity formula for sediment particle, J. Hydraul. Eng., 123(2), 149-152. 
      implicit none
      real(r8), intent(in) :: D_   ! diameter of particle
      real(r8)             :: vsf_       ! sand-sediment erosion rate, [kg/m2/s]
      real(r8) :: R_     ! submerged specific gravity of sediment [-]
      real(r8) :: d_star ! 

      R_ = (densedi-denh2o)/denh2o

      d_star = D_ *((grav*R_/(KVIS**2))**(1._r8/3._r8))
      vsf_ = (KVIS/D_) * (sqrt(25._r8 + 1.2_r8 * d_star*d_star)-5._r8)**1.5_r8

      if(abs(vsf_) < TINYVALUE_s) then
          vsf_ = 0._r8
      endif
      return
    end function CRVSF

  subroutine printTest_sedi(nio)
      ! !DESCRIPTION: output the simulation results into external files
      implicit none
      integer, intent(in) :: nio        ! unit of the file to print

      integer :: IDlist(1:5) = (/1788,2197,1991,3029,1584/)
      integer :: ios,ii                    ! flag of io status


      write(unit=nio,fmt="(25(e20.11))") -TRunoff%etout(IDlist(1),nt_nliq), TRunoff%conc_t(IDlist(1),nt_nsan), TRunoff%wt(IDlist(1),nt_nliq), TRunoff%wt(IDlist(1),nt_nmud), -TRunoff%etout(IDlist(1),nt_nsan), &
                                         -TRunoff%etout(IDlist(2),nt_nliq), TRunoff%conc_t(IDlist(2),nt_nsan), TRunoff%wt(IDlist(2),nt_nliq), TRunoff%wt(IDlist(2),nt_nmud), -TRunoff%etout(IDlist(2),nt_nsan), &
                                         -TRunoff%etout(IDlist(3),nt_nliq), TRunoff%conc_t(IDlist(3),nt_nsan), TRunoff%wt(IDlist(3),nt_nliq), TRunoff%wt(IDlist(3),nt_nmud), -TRunoff%etout(IDlist(3),nt_nsan), &
                                         -TRunoff%etout(IDlist(4),nt_nliq), TRunoff%conc_t(IDlist(4),nt_nsan), TRunoff%wt(IDlist(4),nt_nliq), TRunoff%wt(IDlist(4),nt_nmud), -TRunoff%etout(IDlist(4),nt_nsan), &
                                         -TRunoff%etout(IDlist(5),nt_nliq), TRunoff%conc_t(IDlist(5),nt_nsan), TRunoff%wt(IDlist(5),nt_nliq), TRunoff%wt(IDlist(5),nt_nmud), -TRunoff%etout(IDlist(5),nt_nsan)
  end subroutine printTest_sedi
end MODULE MOSART_sediment_mod