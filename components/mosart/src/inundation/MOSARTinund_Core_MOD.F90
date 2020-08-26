
MODULE MOSARTinund_Core_MOD

!--------------------------------------------------------------------------------------
! DESCRIPTION: Core simulation of MOSART-Inundation.
! 
! HISTORY:
! 2011: Most subroutines/functions were initially developed.
! 2014-2016: The subroutines/functions were revised or created in offline MOSART-Inundation.
! 2017: Integrated with ACME.
! ... ... ...
!
!--------------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_sys_mod, only: shr_sys_abort
  use RunoffMod, only: rtmCTL, Tctl, TUnit, TRunoff, &
      SMatP_upstrm, avsrc_upstrm, avdst_upstrm, SMatP_dnstrm
  use MOSARTinund_PreProcs_MOD, only: con1Em3
  use RtmVar, only: barrier_timers, iulog, inundflag
  use RtmSpmd, only: mpicom_rof, masterproc
  use mct_mod
  
  implicit none
  private
  
  public MOSARTinund_simulate, ManningEq, ChnlFPexchg
                   
  contains
  
  subroutine MOSARTinund_simulate ( )
    ! DESCRIPTION: Simulation of one time step.
    
    ! HISTORY:
    ! 2017: Created for ACME.
    ! ... ... ...
    
    implicit none
    integer :: ier, nr, cnt
    type(mct_aVect) :: avsrc,avdst
    
    ! Hillslope routing computation :
    call inund_hillslopeRouting ( )
    
    ! Subnetwork (tributary channel) routing computation :
    call inund_subnetworkRouting ( )

    ! Store flows from subnetwork to main channel (for budget analysis) :
    TRunoff%erlat_avg = 0._r8
    TRunoff%erlat_avg = TRunoff%erlat_avg - TRunoff%etout            ! Note: outflow is negative for TRunoff%etout .

    if ( Tctl%OPT_inund .eq. 1 ) then
      ! Channel -- floodplain exchange computation :      
      call ChnlFPexchg ( )
    else
      TRunoff%wr_exchg( : ) = TRunoff%wr( :, 1 )
      TRunoff%yr_exchg( : ) = TRunoff%yr( :, 1 )
    end if
    
    if (barrier_timers) then
      call mpi_barrier(mpicom_rof, ier)
    endif

    ! Store channel outflow of previous time step (for budget analysis) :   
    TRunoff%eroup_lagi = 0._r8
    Trunoff%eroup_lagi = Trunoff%eroup_lagi - Trunoff%erout            ! Note: outflow is negative for TRunoff%erout .

    ! Accumulate outflows from upstream channels (for each channel) :
    call AccuUpstrmOutflows ( )

    ! Store outflows from upstream channels (for budget analysis) :
    TRunoff%eroutup_avg = 0._r8
    TRunoff%eroutup_avg = TRunoff%eroutup_avg + TRunoff%eroutUp
    
    ! Channel routing computation :     
    ! Diffusion wave method :
    if ( Tctl%RoutingMethod .eq. 4 ) then               
      
      ! --------------------------------- 
      ! Need code to retrieve values of TRunoff%yr_exchg_dstrm(:) and TRunoff%wr_exchg_dstrm(:) .
      ! --------------------------------- 

      call mct_aVect_init(avsrc,rList='yr:wr',lsize=rtmCTL%lnumr)
      call mct_aVect_init(avdst,rList='yr:wr',lsize=rtmCTL%lnumr)
      call mct_aVect_zero(avsrc)
      call mct_aVect_zero(avdst)
      cnt = 0
      do nr = rtmCTL%begr,rtmCTL%endr
         cnt = cnt + 1
         avsrc%rAttr(1,cnt) = TRunoff%yr_exchg(nr)
         avsrc%rAttr(2,cnt) = TRunoff%wr_exchg(nr)
      enddo

      call mct_sMat_avMult(avsrc, sMatP_dnstrm, avdst)

      cnt = 0
      do nr = rtmCTL%begr,rtmCTL%endr
         cnt = cnt + 1
         TRunoff%yr_exchg_dstrm(nr) = avdst%rAttr(1,cnt)
         TRunoff%wr_exchg_dstrm(nr) = avdst%rAttr(2,cnt)
      enddo

      call mct_aVect_clean(avsrc)
      call mct_aVect_clean(avdst)

      call inund_Routing_DW ( )       
      
    ! Kinematic wave method :
    elseif ( Tctl%RoutingMethod .eq. 1 ) then           
      call inund_Routing_KW ( )       
    end if

    ! Store channel outflow of current time step (for budget analysis) ( actually %flow is same as %eroup_lagf ) :   
    TRunoff%flow = 0._r8
    TRunoff%flow = TRunoff%flow - Trunoff%erout               ! Note: outflow is negative for TRunoff%erout .

    TRunoff%eroup_lagf = 0._r8
    Trunoff%eroup_lagf = Trunoff%eroup_lagf - Trunoff%erout   ! Note: outflow is negative for TRunoff%erout .
    
    ! Transition to the next step :
    call transition2nextStep ( )
  
  end subroutine MOSARTinund_simulate
  
  subroutine inund_hillslopeRouting ( )
    ! DESCRIPTION: Hillslope routing computation.
    
    ! HISTORY:
    ! 2011: Initially developed (H.-Y. Li).
    ! 2014-2016: Revised in offline MOSART-Inundation ( X.Luo ).
    ! 2017: Integrated with ACME.
    ! ... ... ...
    
    implicit none   
    integer :: iu           ! Computation unit index.
    real( r8 ) :: hsV       ! Hillslope flow velocity (m/s).
      
    !$OMP PARALLEL DO PRIVATE(hsV) SCHEDULE(GUIDED)
    do iu = rtmCTL%begr, rtmCTL%endr
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
        
        ! Estimate flow velocity with Manning equation :
        hsV = ManningEq ( TUnit%hslp( iu ), TUnit%nh( iu ), TRunoff%wh( iu, 1 ) )     

        ! Outflow per unit area (outflow = V * A; Negative values means outward flows; Unit: m/s).
        TRunoff%ehout( iu, 1 ) = - hsV * TRunoff%wh( iu, 1 ) * TUnit%Gxr( iu )        

        ! If outflow is too high and water volume becomes negative :
        if ( TRunoff%ehout(iu, 1) .lt. 0._r8 .and. TRunoff%wh(iu, 1) + (TRunoff%qsur(iu, 1) + TRunoff%ehout(iu, 1)) * Tctl%DeltaT .lt. 0._r8 ) then 
          ! All water volume flows away :
          TRunoff%ehout(iu, 1) = - TRunoff%qsur(iu, 1) - TRunoff%wh(iu, 1) / Tctl%DeltaT    
          TRunoff%dwh(iu, 1) = - TRunoff%wh(iu, 1)
        else  
          TRunoff%dwh(iu, 1) = (TRunoff%qsur(iu, 1) + TRunoff%ehout(iu, 1)) * Tctl%DeltaT   ! Added " * Tctl%DeltaT ".
        end if
        
        ! Update water volume over hillslopes :
        TRunoff%wh(iu, 1) = TRunoff%wh(iu, 1) + TRunoff%dwh(iu, 1)    

        ! Total flow rate (surface + subsurface) from hillslopes to subnetwork (tributary channels) within one computation unit  (m^3/s) :
        TRunoff%etin(iu, 1) = ( - TRunoff%ehout(iu, 1) + TRunoff%qsub(iu, 1)) * TUnit%area(iu) * TUnit%frac( iu ) 
      end if  
    end do      
    !$OMP END PARALLEL DO
    
  end subroutine inund_hillslopeRouting
  
  !function CRVRMAN_imp(slp_, n_, rr_) result(v_)
  real( r8 ) function ManningEq( slp_, n_, rr_ )
    ! DESCRIPTION: Estimate flow velocity with Manning equation.
    
    ! HISTORY:
    ! 2011: Initially developed (H.-Y. Li).
    ! 2014-2016: Revised in offline MOSART-Inundation ( X.Luo ).
    ! 2017: Integrated with ACME.
    ! ... ... ...
    
    implicit none
    real( r8 ), intent(in) :: slp_    ! Slope (dimensionless).
    real( r8 ), intent(in) :: n_      ! Manning roughness coefficient ( s * m^(-1/3) ).
    real( r8 ), intent(in) :: rr_     ! Hydraulic radius (m).
    !real( r8 ), intent(out) :: v_    ! Flow velocity (m/s).
    real( r8 ) :: ftemp
    character( len = * ), parameter :: subname = '(ManningEq)'
    
    if ( rr_ .gt. 0._r8 ) then
    
      if ( slp_ .gt. 0._r8 ) then
        ftemp = 2._r8 / 3._r8
        ManningEq = ( rr_ ** ftemp ) * sqrt( slp_ ) / n_     ! Manning equation.
      elseif ( slp_ .lt. 0._r8 ) then
        ftemp = 2._r8 / 3._r8
        ManningEq = - ( rr_ ** ftemp ) * sqrt( - slp_ ) / n_ ! Manning equation.
      else      ! slp_ = 0
        ManningEq = 0._r8
      end if
      
    elseif ( rr_ .eq. 0._r8 ) then
      ManningEq = 0._r8
    else      ! rr_ < 0
      write( iulog, * ) trim( subname ) // ' ERROR: Hydraulic radius is negative !'
      call shr_sys_abort( trim( subname ) // ' ERROR: Hydraulic radius is negative !' )
    end if
    
    return
  end function ManningEq

  subroutine inund_subnetworkRouting ( )
    ! DESCRIPTION: Subnetwork (tributary channel) routing computation.
    
    ! HISTORY:
    ! 2011: Initially developed (H.-Y. Li).
    ! 2014-2016: Revised in offline MOSART-Inundation ( X.Luo ).
    ! 2017: Integrated with ACME.
    ! ... ... ...   
    
    implicit none
    integer :: iu           ! Computation unit index.
    real( r8 ) :: hydrR     ! Hydraulic radius (m).
    
    !$OMP PARALLEL DO PRIVATE(hydrR) SCHEDULE(GUIDED)
    do iu = rtmCTL%begr, rtmCTL%endr
      
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).

        ! If the tributary channel is very short, then no subnetwork routing :
        if ( TUnit%tlen(iu) .lt. Tctl%minL_tribRouting )  then    
          TRunoff%etout(iu, 1) = - TRunoff%etin(iu, 1)    ! Note: outflow is negative.
          TRunoff%dwt(iu, 1) = 0._r8
        else
          ! Calculate hydraulic radius ( hydraulic radius = A / P ) :
          hydrR = ( TUnit%twidth( iu ) * TRunoff%yt( iu, 1 ) ) / ( TUnit%twidth( iu ) + 2._r8 * TRunoff%yt( iu, 1 ) )

          ! Estimate flow velocity with Manning equation :
          TRunoff%vt(iu, 1) = ManningEq (TUnit%tslp(iu), TUnit%nt(iu), hydrR)                 

          ! Calculate flow rate ( = flow velocity * wet cross-sectional area; Note: outflow is negative ) :
          TRunoff%etout(iu, 1) = - TRunoff%vt(iu, 1) * TUnit%twidth( iu ) * TRunoff%yt( iu, 1 )     
          
          ! If outflow rate is too high and water volume becomes negative :
          if ( TRunoff%wt(iu, 1) + (TRunoff%etin(iu, 1) + TRunoff%etout(iu, 1)) * Tctl%DeltaT .lt. 0._r8 ) then 

            ! All water volume flows away :
            TRunoff%etout(iu, 1) = - TRunoff%etin(iu, 1) - TRunoff%wt(iu, 1) / Tctl%DeltaT        

            ! Calculate flow velocity ( = flow rate / wet cross-sectional area ; Note: TRunoff%yt must be non-zero because the initial TRunoff%etout (based on TRunoff%yt) is not zero) :
            TRunoff%vt(iu, 1) = - TRunoff%etout(iu, 1) / TUnit%twidth( iu ) / TRunoff%yt( iu, 1 )
            TRunoff%dwt(iu, 1) = - TRunoff%wt(iu, 1)
          else  
            TRunoff%dwt(iu, 1) = (TRunoff%etin(iu, 1) + TRunoff%etout(iu, 1)) * Tctl%DeltaT   ! Add " * Tctl%DeltaT ".
          end if
        end if
    
        ! Update water volume in subnetwork (tributary channels) :
        TRunoff%wt( iu, 1 ) = TRunoff%wt( iu, 1 ) + TRunoff%dwt( iu, 1 )                

        ! Update water depth :
        if (Tunit%tlen(iu) == 0._r8 .or. Tunit%twidth(iu) == 0._r8) then
           TRunoff%yt(iu,1) = 0._r8
        else
           TRunoff%yt( iu, 1 ) = TRunoff%wt( iu, 1 ) / TUnit%tlen(iu) / TUnit%twidth( iu )     
        endif
      end if
      
    end do
    !$OMP END PARALLEL DO
    
  end subroutine inund_subnetworkRouting

  subroutine ChnlFPexchg ( )
    ! DESCRIPTION: Calculate water exchange between the main channel and floodplains (assuming instantaneous exchange).
    
    ! HISTORY:
    ! 2014-2016: Created and improved in offline MOSART-Inundation ( X.Luo ).
    ! 2017: Integrated with ACME.
    ! ... ... ...

    ! Input:  
    !   + Channel storage and floodplain storage.
    ! Output: 
    !   + New channel storage and floodplain storage after exchange;
    !   + Flooded fraction and flooded area for the floodplains;
    !   + New channel water depth and floodplain max water depth (i.e., elevation difference between final water level (after exchange) and banktop);
    !   + Exchange amount.
    
    implicit none
    integer :: iu
    real( r8 ) :: wr_rcd    ! For recording channel storage before exchange (m^3).    
    real( r8 ) :: w_over    ! = channel storage + floodplain storage - channel storage capacity (m^3).
    integer :: j            ! Index.
    real( r8 ) :: d_s       ! The storage between the current water level and the below water level through the point " j " in the elevation profile (m^3).
    real( r8 ) :: d_e       ! Elevation difference between the current water level and the point " j " in the elevation profile (m).
    real( r8 ) :: hf        ! Floodplain max water depth, namely the elevation difference between the final water level and the banktop (i.e., channel bankfull water level) (m).
    real( r8 ) :: ff_unit   ! Flooded fraction in the computation unit (including channel area) (dimensionless).
    real( r8 ) :: wr_over   ! Channel water storage between the final water level and the banktop ( = channel storage - channel storage capacity ) (m^3).
    character( len = * ), parameter :: subname = '(ChnlFPexchg)'
    
    !$OMP PARALLEL DO PRIVATE(wr_rcd, w_over, j, d_s, d_e, hf, ff_unit, wr_over) SCHEDULE(GUIDED)
    do iu = rtmCTL%begr, rtmCTL%endr
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).

        ! If the channel water level is higher than the floodplain water level, or on the contrary :
        if (TRunoff%yr( iu, 1 ) - TUnit%rdepth(iu) .gt. TRunoff%hf_ini(iu) + con1Em3 .or. &
            (TRunoff%hf_ini(iu) .gt. con1Em3 .and. TRunoff%hf_ini(iu) .gt. TRunoff%yr( iu, 1 ) - TUnit%rdepth(iu) + con1Em3)) then
          wr_rcd = TRunoff%wr( iu, 1 )

          ! Calculate w_over ( = channel storage + floodplain storage - channel storage capacity ) :
          w_over = TRunoff%wr( iu, 1 ) + TRunoff%wf_ini( iu ) - TUnit%wr_bf( iu ) 
        
          ! if ( wr + wf > wr_bf ) :
          if (w_over .gt. 0._r8) then
            ! ---------------------------------  
            ! Calculate: (1) Flooded fraction in the computation unit (including channel area); (2) The difference between final water level and banktop elevation :
            ! --------------------------------- 
            do j = 2, TUnit%npt_eprof3(iu) - 1    ! Note: j starts from 2 .
              if (TUnit%s_eprof3(iu, j) < w_over .and. w_over <= TUnit%s_eprof3(iu, j+1)) then
                ! Water volume above the level of point " j " :
                d_s = w_over - TUnit%s_eprof3(iu, j)

                ! --------------------------------- 
                ! The relationship between d_s and d_e is expressed as a quadratic equation: p3*(d_e)^2 + q3*(d_e) - (d_s) = 0.  
                ! TUnit%p3( :, : ) is the coefficient "p3" of this quadratic equation (unit: m). 
                ! --------------------------------- 
                if (TUnit%p3(iu, j) > 0._r8) then
                  
                  ! Calculate the elevation difference between the current water level and the point " j " ( i.e., the solution of the above quadratic equation ) :
                  d_e = ( - TUnit%q3(iu, j) + sqrt(TUnit%q3(iu, j)**2._r8 + 4._r8*TUnit%p3(iu, j)*d_s)) / (2._r8*TUnit%p3(iu, j))

                  ! Calculate the flooded fraction in the computation unit (including channel area) :
                  ff_unit = d_e * TUnit%alfa3(iu, j) + TUnit%a_eprof3(iu, j)

                ! if ( TUnit%p3(iu, j) = 0 )
                else

                  ! Calculate the elevation difference "d_e" ( = Storage difference / Flooded area ) :
                  !d_e = d_s / TUnit%area(iu) / TUnit%a_eprof3(iu, j)
                  d_e = d_s / ( TUnit%area(iu) * TUnit%frac(iu) * TUnit%a_eprof3(iu, j) )

                  ! Calculate the flooded fraction in the computation unit (including channel area) :
                  ff_unit = TUnit%a_eprof3(iu, j)
                endif

                ! The elevation difference between final water level and banktop :
                hf = TUnit%e_eprof3(iu, j) + d_e
                exit
              endif
            enddo
          
            ! Channel water storage between final water level and banktop :
            wr_over = hf * TUnit%rwidth( iu ) * TUnit%rlen( iu )            

            ! Channel storage after exchange :
            TRunoff%wr_exchg( iu ) = TUnit%wr_bf( iu ) + wr_over        

            ! Floodplain storage after exchange :
            TRunoff%wf_exchg(iu) = w_over - wr_over              
          
            ! Ratio of flooded  area to computation unit area :
            TRunoff%ff_unit(iu) = ff_unit                                                                  
            ! Ratio of flooded floodplain area to computation unit area :
            TRunoff%ff_fp(iu) = ff_unit - TUnit%a_chnl( iu )              

            ! Area of flooded floodplain :
            !TRunoff%fa_fp(iu) = TUnit%area(iu) * TRunoff%ff_fp(iu)      
            TRunoff%fa_fp(iu) = TUnit%area(iu) * TUnit%frac(iu) * TRunoff%ff_fp(iu)

            ! Floodplain max water depth after exchange :
            TRunoff%hf_exchg(iu) = hf                             

            ! Channel water depth after exchange :
            TRunoff%yr_exchg( iu ) = TUnit%rdepth( iu ) + hf            
          
            if ( - con1Em3 < TRunoff%wf_exchg(iu) .and. TRunoff%wf_exchg(iu) < 0._r8) then  
              TRunoff%wf_exchg(iu) = 0._r8            
            elseif ( TRunoff%wf_exchg(iu) <= - con1Em3 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: Floodplain water volume is negative !'
              call shr_sys_abort( trim( subname ) // ' ERROR: Floodplain water volume is negative !' )
            endif   
          
          ! if ( wr + wf <= wr_bf )
          else    
        
            ! The channel water volume after exchange ( all the water remains in the channel ) :
            TRunoff%wr_exchg(iu) = TRunoff%wr(iu,1) + TRunoff%wf_ini(iu)      

            ! The floodplain water volume after exchange is zero :
            TRunoff%wf_exchg(iu) = 0._r8                                  
          
            TRunoff%ff_fp(iu) = 0._r8
            TRunoff%fa_fp(iu) = 0._r8
            TRunoff%hf_exchg(iu) = 0._r8
            TRunoff%ff_unit(iu) = TUnit%a_chnl( iu )

            ! Channel water depth after exchange ( = water volume / channel area ) :
            TRunoff%yr_exchg( iu ) = TRunoff%wr_exchg( iu ) / TUnit%rwidth( iu ) / TUnit%rlen( iu )
          
          endif      ! if ( wr + wf > wr_bf )
        
          ! Channel--floodplain exchange amount (Positive: flow from channel to floodplain; vice versa) (m^3) :
          TRunoff%netchange(iu) = wr_rcd - TRunoff%wr_exchg(iu)
        
        ! ---------------------------------  
        ! No channel--floodplain exchange for two situations: 
        ! (1) Floodplain is inundated: channel water level equals floodplain water level; 
        ! (2) Floodplain is not inundated: channel water level is below banktop.
        ! ---------------------------------  
        else
          TRunoff%wr_exchg( iu ) = TRunoff%wr( iu, 1 )
          TRunoff%yr_exchg( iu ) = TRunoff%yr( iu, 1 )
          TRunoff%wf_exchg( iu ) = TRunoff%wf_ini( iu )
          TRunoff%hf_exchg( iu ) = TRunoff%hf_ini( iu )          
          TRunoff%netchange( iu ) = 0._r8
        end if    ! if ( the channel water level is higher than the floodplain water level, or on the contrary )
      end if      ! if ( TUnit%mask( iu ) .gt. 0 )
    end do
    !$OMP END PARALLEL DO
    
  end subroutine ChnlFPexchg
  
  subroutine AccuUpstrmOutflows ( )
    ! DESCRIPTION: Accumulate outflows from upstream channels (for each channel).
    
    implicit none
    integer :: iu, k, cnt
    integer :: unitUp       ! ID of upstream channel.
    
    TRunoff%eroutUp = 0._r8
    
#ifdef NO_MCT
    do iu = rtmCTL%begr, rtmCTL%endr
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
      
        do k=1, TUnit%nUp(iu)
          unitUp = TUnit%iUp(iu, k)
          TRunoff%eroutUp(iu, 1) = TRunoff%eroutUp(iu, 1) + TRunoff%erout(unitUp, 1)    ! Note: outflow is negative for TRunoff%erout .
        end do
        
      end if  
    end do
#else
    !--- copy erout into avsrc_upstrm ---
    call mct_avect_zero( avsrc_upstrm )
    cnt = 0
    do iu = rtmCTL%begr, rtmCTL%endr
      cnt = cnt + 1     
      avsrc_upstrm%rAttr(1, cnt) = TRunoff%erout(iu, 1)
    enddo
    
    call mct_avect_zero( avdst_upstrm )

    call mct_sMat_avMult(avsrc_upstrm, sMatP_upstrm, avdst_upstrm)

    !--- add mapped eroutUp to TRunoff ---
    cnt = 0
    do iu = rtmCTL%begr, rtmCTL%endr
      cnt = cnt + 1     
      TRunoff%eroutUp(iu, 1) = avdst_upstrm%rAttr(1, cnt)
    enddo
#endif    

  end subroutine AccuUpstrmOutflows
  
  subroutine inund_Routing_DW ( )
    ! DESCRIPTION: Channel routing computation with diffusion wave method.
    
    ! HISTORY:
    ! 2011: Initially developed (H.-Y. Li).
    ! 2014-2016: Revised in offline MOSART-Inundation ( X.Luo ).
    ! 2017: Integrated with ACME.
    ! ... ... ...
    
    implicit none
    integer :: iu
    real( r8 ) :: surfaceSlope                  ! Water surface slope ( from the current channel to the downstream channel ) (dimensionless).
    real( r8 ) :: y_c, len_c, slp_c             ! Water depth (m), length (m) and bed slope (dimensionless) of the current channel.
    real( r8 ) :: y_down, len_down, slp_down    ! Water depth (m), length (m) and bed slope (dimensionless) of the downstream channel.
    real( r8 ) :: hydrR                         ! Hydraulic radius (m).   
    real( r8 ) :: wv_gwl                        ! Water volume from glacier, wetlands or lakes (m^3). 
    character( len = * ), parameter :: subname = '(inund_Routing_DW)'
    
    !$OMP PARALLEL DO PRIVATE(surfaceSlope, y_c, len_c, slp_c, y_down, len_down, slp_down, hydrR, wv_gwl) SCHEDULE(GUIDED)
    do iu = rtmCTL%begr, rtmCTL%endr      
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).

        ! ---------------------------------   
        ! Calculate water surface slope ( from the current channel to the downstream channel ) :
        ! ---------------------------------  

        ! The water surface slope is used to estimate flow velocity :
        if ( Tctl%OPT_trueDW .eq. 1 ) then

          ! If this channel is NOT at basin outlet :
          !if ( TUnit%mask( iu ) .eq. 1 ) then     
          if ( rtmCTL%mask( iu ) .eq. 1 ) then
            y_c = TRunoff%yr_exchg( iu )
            len_c = TUnit%rlen( iu )
            slp_c = TUnit%rslp( iu )
            y_down = TRunoff%yr_exchg_dstrm( iu )
            len_down = TUnit%rlen_dstrm( iu )
            slp_down = TUnit%rslp_dstrm( iu )

          ! If this channel is at basin outlet (downstream is ocean) :
          !elseif ( TUnit%mask( iu ) .eq. 2 ) then
          elseif ( rtmCTL%mask( iu ) .eq. 3 ) then
            y_c = TRunoff%yr_exchg( iu )
            len_c = TUnit%rlen( iu )
            slp_c = TUnit%rslp( iu )

            ! ---------------------------------  
            ! Water depth equals channel depth at the basin outlet while assuming :
            ! (1) The sea level equals the channel banktop, and is fixed;
            ! (2) The river stage equals the sea level. 
            ! ---------------------------------  
            y_down = TUnit%rdepth( iu )

            len_down = 0._r8
            slp_down = 0._r8    
          end if
      
          ! Calculate water surface slope ( from current-channel surface mid-point to downstream-channel surface mid-point ) :
          surfaceSlope = (len_down * slp_down + len_c * slp_c + 2._r8 * y_c - 2._r8 * y_down) / (len_c + len_down)

        ! Riverbed slope is used as the surrogate for water surface slope ( this is a temporary treatment before the downstream-channel information can be retrieved ) :
        elseif ( Tctl%OPT_trueDW .eq. 2 ) then
        
          surfaceSlope = TUnit%rslp( iu )

        endif
      
        ! ----------------------------------   
        ! Calculate flow velocity, streamflow and water volume change :
        ! ---------------------------------  

        if ( surfaceSlope .gt. 0._r8 ) then

          ! If water volume becomes negative because upward flow (from current channel to upstream channels) is too large :
          if ( TRunoff%wr_exchg( iu ) + ( - TRunoff%etout( iu, 1 ) - TRunoff%eroutUp( iu, 1 ) ) * Tctl%DeltaT .lt. 0._r8 ) then
            TRunoff%vr( iu, 1 ) = 0._r8
            TRunoff%erout( iu, 1 ) = 0._r8
            TRunoff%dwr( iu, 1 ) = - TRunoff%wr_exchg( iu )       

          else
            ! Hydraulic radius ( = A / P ) :
            hydrR = ( TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu ) ) / ( TUnit%rwidth( iu ) + 2._r8 * TRunoff%yr_exchg( iu ) )

            ! Estimate flow velocity with Manning equation :
            TRunoff%vr( iu, 1 ) = ManningEq ( surfaceSlope, TUnit%nr( iu ), hydrR )
          
            ! Calculate outflow rate ( Q = V * A ; Note: outflow is negative) :
            TRunoff%erout( iu, 1 ) = - TRunoff%vr( iu, 1 ) * TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu )
          
            ! Water volume change :
            TRunoff%dwr( iu, 1 ) = ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) + TRunoff%erout( iu, 1 ) ) * Tctl%DeltaT
          
            ! If water volume becomes negative because downward flow is too large :
            if ( TRunoff%erout( iu, 1 ) .lt. 0._r8 .and. TRunoff%wr_exchg( iu ) + TRunoff%dwr( iu, 1 ) .lt. 0._r8 ) then
            
              ! All the water flows away :
              TRunoff%erout( iu, 1 ) = - ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) + TRunoff%wr_exchg( iu ) / Tctl%DeltaT )

              if ( TRunoff%yr_exchg( iu ) .gt. 0._r8 ) then
                TRunoff%vr( iu, 1 ) = - TRunoff%erout( iu, 1 ) / TUnit%rwidth( iu ) / TRunoff%yr_exchg( iu )
              else
                write( iulog, * ) trim( subname ) // ' ERROR: Downward flow is too large & water depth is zero !'
                call shr_sys_abort( trim( subname ) // ' ERROR: Downward flow is too large & water depth is zero !' )
              end if
              TRunoff%dwr( iu, 1 ) = - TRunoff%wr_exchg( iu )
            end if
          end if
        
        elseif ( surfaceSlope .lt. 0._r8 ) then

          ! Hydraulic radius ( = A / P ) :
          hydrR = ( TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu ) ) / ( TUnit%rwidth( iu ) + 2._r8 * TRunoff%yr_exchg( iu ) )

          ! Estimate flow velocity with Manning equation ( flow velocity is negative, i.e., upward flow ) :
          TRunoff%vr( iu, 1 ) = ManningEq ( surfaceSlope, TUnit%nr( iu ), hydrR )
        
          ! Calculate outflow rate ( Q = V * A ; Flow rate is positive, i.e., upward flow ) :
          TRunoff%erout( iu, 1 ) = - TRunoff%vr( iu, 1 ) * TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu )
        
          ! If this channel is NOT at basin outlet :
          !if ( TUnit%mask( iu ) .eq. 1 ) then
          if ( rtmCTL%mask( iu ) .eq. 1 ) then
          
            ! If upward flow (from downstream channel to current channel) is too large, so that water volume becomes negative in downstream channel :
            if ( TRunoff%wr_exchg_dstrm( iu ) - TRunoff%erout( iu, 1 ) * Tctl%DeltaT .lt. 0._r8 ) then
            
              TRunoff%erout( iu, 1 ) = TRunoff%wr_exchg_dstrm( iu ) / Tctl%DeltaT
              if ( TRunoff%yr_exchg( iu ) .gt. 0._r8 ) then
                TRunoff%vr( iu, 1 ) = - TRunoff%erout( iu, 1 ) / TUnit%rwidth( iu ) / TRunoff%yr_exchg( iu )
              else
                write( iulog, * ) trim( subname ) // ' ERROR: Upward flow is too large & water depth is zero !'
                call shr_sys_abort( trim( subname ) // ' ERROR: Upward flow is too large & water depth is zero !' )
              end if
            
              !TRunoff%wr_exchg_dstrm( iu ) = 0._r8     ! Because all downstream-channel water flows to the current channel.
            end if
          end if
        
          ! Water volume change :
          TRunoff%dwr( iu, 1 ) = ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) + TRunoff%erout( iu, 1 ) ) * Tctl%DeltaT
        
        ! If water surface slope (surfaceSlope) is zero :
        else
          TRunoff%vr( iu, 1 ) = 0._r8
          TRunoff%erout( iu, 1 ) = 0._r8
          TRunoff%dwr( iu, 1 ) = ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) ) * Tctl%DeltaT
        end if
      
        ! Channel water volume after routing computation :
        TRunoff%wr_rtg( iu ) = TRunoff%wr_exchg( iu ) + TRunoff%dwr( iu, 1 )

        ! Account for TRunoff%qgwl(...) ( runoff from glacier, wetlands or lakes; unit: m/s ) :
        wv_gwl = TRunoff%qgwl( iu, 1 ) * TUnit%area( iu ) * TUnit%frac( iu ) * Tctl%DeltaT

        if ( TRunoff%wr_rtg( iu ) + wv_gwl .lt. 0._r8 ) then
          write( iulog, * ) trim( subname ) // ' WARNING: TRunoff%qgwl(...) is so negative that the channel water volume becomes negative ! Channel water volume is forced to be zero.'
          write( iulog, * ) trim( subname ) // ' WARNING: at the computation unit ', iu

          TRunoff%wr_rtg( iu ) = 0._r8
        end if  

        ! Channel water depth after routing computation :
        TRunoff%yr_rtg( iu ) = TRunoff%wr_rtg( iu ) / TUnit%rwidth( iu ) / TUnit%rlen( iu )
    
      end if    ! end of if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 )
    end do
    !$OMP END PARALLEL DO
    
  end subroutine inund_Routing_DW
  
  subroutine inund_Routing_KW ( )
    ! DESCRIPTION: Channel routing computation with kinematic wave method.
    
    ! HISTORY:
    ! 2011: Initially developed (H.-Y. Li).
    ! 2014-2016: Revised in offline MOSART-Inundation ( X.Luo ).
    ! 2017: Integrated with ACME.
    ! ... ... ...
    
    implicit none
    integer :: iu
    real( r8 ) :: hydrR     ! Hydraulic radius (m).
    real( r8 ) :: wv_gwl    ! Water volume from glacier, wetlands or lakes (m^3).
    character( len = * ), parameter :: subname = '(inund_Routing_KW)'
    
    !$OMP PARALLEL DO PRIVATE(hydrR, wv_gwl) SCHEDULE(GUIDED)
    do iu = rtmCTL%begr, rtmCTL%endr
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).      

        ! ---------------------------------  
        ! Calculate flow velocity, streamflow and water volume change :
        ! ---------------------------------  

        ! Hydraulic radius ( = A / P ) :
        hydrR = ( TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu ) ) / ( TUnit%rwidth( iu ) + 2._r8 * TRunoff%yr_exchg( iu ) )

        ! Estimate flow velocity with Manning equation :
        TRunoff%vr( iu, 1 ) = ManningEq ( TUnit%rslp( iu ), TUnit%nr( iu ), hydrR )     
      
        ! Calculate outflow rate ( Q = V * A ; Note: outflow is negative ) :
        TRunoff%erout( iu, 1 ) = - TRunoff%vr( iu, 1 ) * TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu )      
      
        ! Water volume change :
        TRunoff%dwr( iu, 1 ) = ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) + TRunoff%erout( iu, 1 ) ) * Tctl%DeltaT   
      
        ! If water volume becomes negative because downward flow is too large :
        if ( TRunoff%erout( iu, 1 ) .lt. 0._r8 .and. TRunoff%wr_exchg( iu ) + TRunoff%dwr( iu, 1 ) .lt. 0._r8 ) then    
          TRunoff%erout( iu, 1 ) = - ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) + TRunoff%wr_exchg( iu ) / Tctl%DeltaT )
          if ( TRunoff%yr_exchg( iu ) .gt. 0._r8 ) then
            TRunoff%vr( iu, 1 ) = - TRunoff%erout( iu, 1 ) / TUnit%rwidth( iu ) / TRunoff%yr_exchg( iu )
          else
            write( iulog, * ) trim( subname ) // ' ERROR: Downward flow is too large & water depth is zero !'
            call shr_sys_abort( trim( subname ) // ' ERROR: Downward flow is too large & water depth is zero !' )
          end if
          TRunoff%dwr( iu, 1 ) = - TRunoff%wr_exchg( iu )
        end if
      
        ! Channel water volume after routing computation :
        TRunoff%wr_rtg( iu ) = TRunoff%wr_exchg( iu ) + TRunoff%dwr( iu, 1 )    

        ! Account for TRunoff%qgwl(...) ( runoff from glacier, wetlands or lakes; unit: m/s ) :
        wv_gwl = TRunoff%qgwl( iu, 1 ) * TUnit%area( iu ) * TUnit%frac( iu ) * Tctl%DeltaT

        if ( TRunoff%wr_rtg( iu ) + wv_gwl .lt. 0._r8 ) then
          write( iulog, * ) trim( subname ) // ' WARNING: TRunoff%qgwl(...) is so negative that the channel water volume becomes negative ! Channel water volume is forced to be zero.'
          write( iulog, * ) trim( subname ) // ' WARNING: at the computation unit ', iu

          TRunoff%wr_rtg( iu ) = 0._r8
        end if  

        ! Channel water depth after routing computation :
        TRunoff%yr_rtg( iu ) = TRunoff%wr_rtg( iu ) / TUnit%rwidth( iu ) / TUnit%rlen( iu )   
    
      end if    ! end of if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 )
    end do
    !$OMP END PARALLEL DO
    
  end subroutine inund_Routing_KW
  
  subroutine transition2nextStep ( )
    ! DESCRIPTION: Transition to the next step ( transfer values to variables of the next step ).
    implicit none
    
    ! Channel water volume :
    TRunoff%wr( :, 1 ) = TRunoff%wr_rtg( : )      

    ! Channel water depth :
    TRunoff%yr( :, 1 ) = TRunoff%yr_rtg( : )      

    ! Channel streamflow :
    !TRunoff%chnlQ_prev = TRunoff%chnlQ   
 
    if ( Tctl%OPT_inund .eq. 1 ) then
      ! Floodplain water volume :
      TRunoff%wf_ini = TRunoff%wf_exchg     

      ! Floodplain max water depth :
      TRunoff%hf_ini = TRunoff%hf_exchg     
      ! Floodplain area fraction (not including channel)
      TRunoff%ff_ini = TRunoff%ff_fp
      ! Flooded area fraction (including channel):
      TRunoff%ffunit_ini = TRunoff%ff_unit
    end if    
  
  end subroutine transition2nextStep

!#endif
  
end MODULE MOSARTinund_Core_MOD
