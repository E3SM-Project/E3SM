
MODULE MOSARTinund_Core_MOD

!------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
! DESCRIPTION: Core simulation of MOSART-Inundation.
! 
! HISTORY:
! 2011: Most subroutines/functions were initially developed.
! 2014-2016: The subroutines/functions were revised or created in offline MOSART-Inundation.
! 2017: Integrated with ACME.
! ... ... ...
!
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	use shr_kind_mod_inund, only: r8 => shr_kind_r8
	use RunoffMod, only: rtmCTL, Tctl, TUnit, TRunoff, SMatP_eroutUp, avsrc_eroutUp, avdst_eroutUp
	use MOSARTinund_PreProcs_MOD, only: con1Em3
	use RtmVar, only: barrier_timers
	use RtmSpmd, only: mpicom_rof
	
	implicit none
	
	public MOSARTinund_simulate
	private inund_hillslopeRouting, ManningEq, inund_subnetworkRouting, ChnlFPexchg, inund_Routing_DW, inund_Routing_KW, transition2nextStep
								   
	contains
	
	subroutine MOSARTinund_simulate ( )
		! DESCRIPTION: Simulation of one time step.
		
		! HISTORY:
		! 2017: Created for ACME.
		! ... ... ...
		
		implicit none
		integer :: ier, iu, k, cnt
		integer :: unitUp			! ID of upstream channel.
		
		! Hillslope routing computation :
		call inund_hillslopeRouting ( )
		
		! Subnetwork (tributary channel) routing computation :
		call inund_subnetworkRouting ( )
		
		if ( Tctl%OPT_inund .eq. 1 ) then
			! Channel -- floodplain exchange computation :			
			call ChnlFPexchg ( )
		else
			TRunoff%wr_exchg( : ) = TRunoff%wr( :, 1 )
			TRunoff%yr_exchg( : ) = TRunoff%yr( :, 1 )
		end if
		
		! Accumulate outflows from upstream channels (for each channel) (i.e., TRunoff%eroutUp(:, :) ) :
		
		if (barrier_timers) then
			call mpi_barrier(mpicom_rof, ier)
		endif

		TRunoff%eroup_lagi = 0._r8
		
		! Store channel outflow of previous time step (for budget analysis) :		
		Trunoff%eroup_lagi = Trunoff%eroup_lagi - Trunoff%erout

		TRunoff%eroutUp = 0._r8
#ifdef NO_MCT
		do iu = rtmCTL%begr, rtmCTL%endr
			if ( TUnit%mask( iu ) .gt. 0 ) then
			
				do k=1, TUnit%nUp(iu)
					unitUp = TUnit%iUp(iu, k)
					TRunoff%eroutUp(iu, 1) = TRunoff%eroutUp(iu, 1) + TRunoff%erout(unitUp, 1)
				end do
				
			end if	
		end do
#else
		!--- copy erout into avsrc_eroutUp ---
		call mct_avect_zero( avsrc_eroutUp )
		cnt = 0
		do iu = rtmCTL%begr, rtmCTL%endr
			cnt = cnt + 1			
			avsrc_eroutUp%rAttr(1, cnt) = TRunoff%erout(iu, 1)
		enddo
		
		call mct_avect_zero( avdst_eroutUp )

		call mct_sMat_avMult(avsrc_eroutUp, sMatP_eroutUp, avdst_eroutUp)

		!--- add mapped eroutUp to TRunoff ---
		cnt = 0
		do iu = rtmCTL%begr, rtmCTL%endr
			cnt = cnt + 1			
			TRunoff%eroutUp(iu, 1) = avdst_eroutUp%rAttr(1, cnt)
		enddo
#endif		
		
		! Channel routing computation :			
		if ( Tctl%RoutingMethod .eq. 4 ) then				! Diffusion wave method.				
			
			! --------------------------------- 
			! Need code to retrieve values of TRunoff%yr_exchg_dstrm(:) and TRunoff%wr_exchg_dstrm(:) .
			! --------------------------------- 
		
			call inund_Routing_DW ( )				
			
		elseif ( Tctl%RoutingMethod .eq. 1 ) then		! Kinematic wave method.				
			call inund_Routing_KW ( )				
		end if
		
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
		integer :: iu								! Computation unit index.
		real( r8 ) :: hsV						! Hillslope flow velocity (m/s).
			
		!$OMP PARALLEL FIRSTPRIVATE( iu, hsV )
		!$OMP DO SCHEDULE( GUIDED )
		do iu = rtmCTL%begr, rtmCTL%endr
			if ( TUnit%mask( iu ) .gt. 0 ) then
				!TRunoff%ehout(iu) = -CREHT_imp(TUnit%hslp_input(iu), TUnit%nh_calc(iu), TUnit%Gxr_input(iu), TRunoff%wh(iu))		! (m/s) ('Negative' means flow out of hillslope.) ; Replace %yh with %wh.
				hsV = ManningEq ( TUnit%hslp( iu ), TUnit%nh( iu ), TRunoff%wh( iu, 1 ) )			! Estimate flow velocity with Manning equation.
				TRunoff%ehout( iu, 1 ) = - hsV * TRunoff%wh( iu, 1 ) * TUnit%Gxr( iu )				! Outflow per unit area (outflow = V * A; Negative values means outward flows) (m/s).
				
				if ( TRunoff%ehout(iu, 1) .lt. 0._r8 .and. TRunoff%wh(iu, 1) + (TRunoff%qsur(iu, 1) + TRunoff%ehout(iu, 1)) * Tctl%inund_dt .lt. 0._r8 ) then	! Outflow is too high and water volume becomes negative.
					TRunoff%ehout(iu, 1) = - TRunoff%qsur(iu, 1) - TRunoff%wh(iu, 1) / Tctl%inund_dt  	! All water volume flows away.
					TRunoff%dwh(iu, 1) = - TRunoff%wh(iu, 1)
				else	
					TRunoff%dwh(iu, 1) = (TRunoff%qsur(iu, 1) + TRunoff%ehout(iu, 1)) * Tctl%inund_dt   ! Added " * Tctl%inund_dt ".
				end if
				
				TRunoff%wh(iu, 1) = TRunoff%wh(iu, 1) + TRunoff%dwh(iu, 1)		! Update water volume.
				TRunoff%etin(iu, 1) = ( - TRunoff%ehout(iu, 1) + TRunoff%qsub(iu, 1)) * TUnit%area(iu) * TUnit%frac( iu )	! Total flow rate (surface + subsurface) from hillslopes to subnetwork (tributary channels) within one computation unit  (m^3/s) .
			end if	
		end do			
		!$OMP END DO
		!$OMP END PARALLEL	
		
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
		real( r8 ), intent(in) :: slp_ 		! Slope (dimensionless).
		real( r8 ), intent(in) :: n_ 			! Manning roughness coefficient ( s * m^(-1/3) ).
		real( r8 ), intent(in) :: rr_ 		! Hydraulic radius (m).
		!real( r8 ), intent(out) :: v_     	! Flow velocity (m/s).
		real( r8 ) :: ftemp
		
		if ( rr_ .gt. 0._r8 ) then
		
			if ( slp_ .gt. 0._r8 ) then
				ftemp = 2._r8 / 3._r8
				ManningEq = ( rr_ ** ftemp ) * sqrt( slp_ ) / n_			! Manning equation.
			elseif ( slp_ .lt. 0._r8 ) then
				ftemp = 2._r8 / 3._r8
				ManningEq = - ( rr_ ** ftemp ) * sqrt( - slp_ ) / n_		! Manning equation.
			else			! slp_ = 0
				ManningEq = 0._r8
			end if
			
		elseif ( rr_ .eq. 0._r8 ) then
			ManningEq = 0._r8
		else			! rr_ < 0
			write ( *, * ) "Function ManningEq ( ): Hydraulic radius is negative !"
			STOP
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
		integer :: iu									! Computation unit index.
		real( r8 ) :: hydrR						! Hydraulic radius (m).
		
		!$OMP PARALLEL FIRSTPRIVATE ( iu, hydrR )
		!$OMP DO SCHEDULE ( GUIDED )
		do iu = rtmCTL%begr, rtmCTL%endr
			
			if ( TUnit%mask( iu ) .gt. 0 ) then
				if ( TUnit%tlen(iu) .lt. Tctl%minL_tribRouting )	then		! If the tributary channel is very short, then no subnetwork routing.
					TRunoff%etout(iu, 1) = - TRunoff%etin(iu, 1)			! Note: outflow is negative.
					TRunoff%dwt(iu) = 0._r8
				else
					hydrR = ( TUnit%twidth( iu ) * TRunoff%yt( iu, 1 ) ) / ( TUnit%twidth( iu ) + 2._r8 * TRunoff%yt( iu, 1 ) )		! hydrR = A / P.
					TRunoff%vt(iu, 1) = ManningEq (TUnit%tslp(iu), TUnit%nt(iu), hydrR)									! Estimate flow velocity with Manning equation.
					TRunoff%etout(iu, 1) = - TRunoff%vt(iu, 1) * TUnit%twidth( iu )	* TRunoff%yt( iu, 1 )			! Flow rate = flow velocity * wet cross-sectional area.
					
					if ( TRunoff%wt(iu, 1) + (TRunoff%etin(iu, 1) + TRunoff%etout(iu, 1)) * Tctl%inund_dt .lt. 0._r8 ) then	! Outflow rate is too high and water volume becomes negative.
						TRunoff%etout(iu, 1) = - TRunoff%etin(iu, 1) - TRunoff%wt(iu, 1) / Tctl%inund_dt  			! All water volume flows away.
						TRunoff%vt(iu, 1) = - TRunoff%etout(iu, 1) / TUnit%twidth( iu )	/ TRunoff%yt( iu, 1 )		! Flow velocity = flow rate / wet cross-sectional area ( Note: TRunoff%yt cannot be zero because TRunoff%etout is not zero ).
						TRunoff%dwt(iu, 1) = - TRunoff%wt(iu, 1)
					else	
						TRunoff%dwt(iu, 1) = (TRunoff%etin(iu, 1) + TRunoff%etout(iu, 1)) * Tctl%inund_dt   ! Add " * Tctl%inund_dt ".
					end if
				end if
		
				TRunoff%wt( iu, 1 ) = TRunoff%wt( iu, 1 ) + TRunoff%dwt( iu, 1 )								! Update water volume.
				TRunoff%yt( iu, 1 ) = TRunoff%wt( iu, 1 ) / TUnit%tlen(iu) / TUnit%twidth( iu )			! Update water depth.
			end if
			
		end do
		!$OMP END DO
		!$OMP END PARALLEL	
		
	end subroutine inund_subnetworkRouting

	subroutine ChnlFPexchg ( )
		! DESCRIPTION: Calculate water exchange between the main channel and floodplains (assuming instantaneous exchange).
		
		! HISTORY:
		! 2014-2016: Created and improved in offline MOSART-Inundation ( X.Luo ).
		! 2017: Integrated with ACME.
		! ... ... ...

		! Input:	
		!		+ Channel storage and floodplain storage.
		! Output:	
		!		+ New channel storage and floodplain storage after exchange;
		!		+ Flooded fraction and flooded area for the floodplains;
		!		+ New channel water depth and floodplain max water depth (i.e., elevation difference between final water level (after exchange) and banktop);
		!		+ Exchange amount.
		
		implicit none
		integer :: iu
		real( r8 ) :: wr_rcd		! For recording channel storage before exchange (m^3).		
		real( r8 ) :: w_over		! = channel storage + floodplain storage - channel storage capacity (m^3).
		integer :: j					! Index.
		real( r8 ) :: d_s			! The storage between the current water level and the below water level through the point " j " in the elevation profile (m^3).
		real( r8 ) :: d_e			! Elevation difference between the current water level and the point " j " in the elevation profile (m).
		real( r8 ) :: hf				! Floodplain max water depth, namely the elevation difference between the final water level and the banktop (i.e., channel bankfull water level) (m).
		real( r8 ) :: ff_unit		! Flooded fraction in the computation unit (including channel area) (dimensionless).
		real( r8 ) :: wr_over	! Channel water storage between the final water level and the banktop ( = channel storage - channel storage capacity ) (m^3).
		
		!$OMP PARALLEL FIRSTPRIVATE( iu, wr_rcd, w_over, j, d_s, d_e, hf, ff_unit, wr_over )
		!$OMP DO SCHEDULE( GUIDED )
		do iu = rtmCTL%begr, rtmCTL%endr
			if ( TUnit%mask( iu ) .gt. 0 ) then
				if (TRunoff%yr( iu, 1 ) - TUnit%rdepth(iu) .gt. TRunoff%hf_ini(iu) + con1Em3 .or. &
						(TRunoff%hf_ini(iu) .gt. con1Em3 .and. TRunoff%hf_ini(iu) .gt. TRunoff%yr( iu, 1 ) - TUnit%rdepth(iu) + con1Em3)) then		! The channel water level is higher than the floodplain water level, or on the contrary.
					wr_rcd = TRunoff%wr( iu, 1 )
					w_over = TRunoff%wr( iu, 1 ) + TRunoff%wf_ini( iu ) - TUnit%wr_bf( iu )	! = channel storage + floodplain storage - channel storage capacity
				
					if (w_over .gt. 0._r8) then 		! wr + wf > wr_bf
						! Calculate the difference between final water level and banktop elevation :
						do j = 2, TUnit%npt_eprof3(iu) - 1		! NOTE: j starts from 2 .
							if (TUnit%s_eprof3(iu, j) < w_over .and. w_over <= TUnit%s_eprof3(iu, j+1)) then
								d_s = w_over - TUnit%s_eprof3(iu, j)
								if (TUnit%p3(iu, j) > 0._r8) then				! '%p3' is a coefficient (unit: m) of the following quadratic equation.
									! The relationship between d_s and d_e is expressed as a quadratic equation: p3*(d_e)^2 + q3*(d_e) - (d_s) = 0. 
									d_e = ( - TUnit%q3(iu, j) + sqrt(TUnit%q3(iu, j)**2._r8 + 4._r8*TUnit%p3(iu, j)*d_s)) / (2._r8*TUnit%p3(iu, j))	! Elevation difference is the solution of the quadratic equation.
									ff_unit = d_e * TUnit%alfa3(iu, j) + TUnit%a_eprof3(iu, j)	! Flooded fraction in the computation unit (including channel area).
								else		! p3 = 0
									d_e = d_s / TUnit%area(iu) / TUnit%a_eprof3(iu, j)	! Elevation difference = Storage difference / Flooded area
									ff_unit = TUnit%a_eprof3(iu, j)
								endif
								hf = TUnit%e_eprof3(iu, j) + d_e		! Elevation difference between final water level and banktop.
								exit
							endif
						enddo
					
						wr_over = hf * TUnit%rwidth( iu ) * TUnit%rlen( iu )						! Channel water storage between final water level and banktop.
						TRunoff%wr_exchg( iu ) = TUnit%wr_bf( iu ) + wr_over				! Channel storage after exchange.
						TRunoff%wf_exchg(iu) = w_over - wr_over								! Floodplain storage after exchange.
					
						TRunoff%ff_fp(iu) = ff_unit - TUnit%a_chnl( iu )							! = area of flooded floodplain / unit area .
						TRunoff%fa_fp(iu) = TUnit%area(iu) * TRunoff%ff_fp(iu)			! Area of flooded floodplain .
						TRunoff%hf_exchg(iu) = hf															! Floodplain max water depth after exchange .
						TRunoff%yr_exchg( iu ) = TUnit%rdepth( iu ) + hf						! Channel water depth after exchange.
					
						if ( - con1Em3 < TRunoff%wf_exchg(iu) .and. TRunoff%wf_exchg(iu) < 0._r8) then	
							TRunoff%wf_exchg(iu) = 0._r8						
						elseif ( TRunoff%wf_exchg(iu) <= - con1Em3 ) then
							write (*, *) "Subroutine ChnlFPexchg ( ) : Floodplain water volume is NEGATIVE !"
							STOP
						endif		
					
					else		! wr + wf <= wr_bf
				
						TRunoff%wr_exchg(iu) = TRunoff%wr(iu,1) + TRunoff%wf_ini(iu)			! All the water remains in the channel.
						TRunoff%wf_exchg(iu) = 0._r8																	! Floodplain is not inundated.
					
						TRunoff%ff_fp(iu) = 0._r8
						TRunoff%fa_fp(iu) = 0._r8
						TRunoff%hf_exchg(iu) = 0._r8
						TRunoff%yr_exchg( iu ) = TRunoff%wr_exchg( iu ) / TUnit%rwidth( iu ) / TUnit%rlen( iu )
					
					endif	! if ( wr + wf > wr_bf )
				
					TRunoff%se_rf(iu) = wr_rcd - TRunoff%wr_exchg(iu)	! Channel -- floodplain exchange amount (Positive: flow from channel to floodplain; vice versa) (m^3). 
				
				else			! No channel--floodplain exchange for two situations: (1) Floodplain is inundated: channel water level equals floodplain water level; (2) Floodplain is not inundated: channel water level is below banktop.
					TRunoff%wr_exchg( iu ) = TRunoff%wr( iu, 1 )
					TRunoff%wf_exchg( iu ) = TRunoff%wf_ini( iu )
					TRunoff%hf_exchg( iu ) = TRunoff%hf_ini( iu )
					TRunoff%yr_exchg( iu ) = TRunoff%yr( iu, 1 )
					TRunoff%se_rf( iu ) = 0._r8
				end if		! if ( the channel water level is higher than the floodplain water level, or on the contrary )
			end if 		! if ( TUnit%mask( iu ) .gt. 0 )
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		
	end subroutine ChnlFPexchg
	
	subroutine inund_Routing_DW ( )
		! DESCRIPTION: Channel routing computation with diffusion wave method.
		
		! HISTORY:
		! 2011: Initially developed (H.-Y. Li).
		! 2014-2016: Revised in offline MOSART-Inundation ( X.Luo ).
		! 2017: Integrated with ACME.
		! ... ... ...
		
		implicit none
		integer :: iu
		integer :: k
		real( r8 ) :: surfaceSlope			! Water surface slope ( from the current channel to the downstream channel ) (dimensionless).
		real( r8 ) :: y_c, len_c, slp_c		! Water depth (m), length (m) and bed slope (dimensionless) of the current channel.
		real( r8 ) :: y_down, len_down, slp_down		! Water depth (m), length (m) and bed slope (dimensionless) of the downstream channel.
		real( r8 ) :: hydrR						! Hydraulic radius (m).		
		
		!$OMP PARALLEL FIRSTPRIVATE( iu, k, surfaceSlope, y_c, len_c, slp_c, y_down, len_down, slp_down, hydrR )
		!$OMP DO SCHEDULE( GUIDED )
		do iu = rtmCTL%begr, rtmCTL%endr			
			if ( TUnit%mask( iu ) .gt. 0 ) then
				! Calculate water surface slope ( from the current channel to the downstream channel ) :
				if ( TUnit%mask( iu ) .eq. 1 ) then			! This channel is NOT at basin outlet.
					y_c = TRunoff%yr_exchg( iu )
					len_c = TUnit%rlen( iu )
					slp_c = TUnit%rslp( iu )
					y_down = TRunoff%yr_exchg_dstrm( iu )
					len_down = TUnit%rlen_dstrm( iu )
					slp_down = TUnit%rslp_dstrm( iu )
				elseif ( TUnit%mask( iu ) .eq. 2 ) then			! This channel is at basin outlet (downstream is ocean).
					y_c = TRunoff%yr_exchg( iu )
					len_c = TUnit%rlen( iu )
					slp_c = TUnit%rslp( iu )
					y_down = TUnit%rdepth( iu )		! Water depth equals channel depth at the basin outlet while assuming : (1) The sea level equals the channel banktop, and is fixed; (2) The river stage equals the sea level. 
					len_down = 0._r8
					slp_down = 0._r8		
				end if
			
				surfaceSlope = (len_down * slp_down + len_c * slp_c + 2 * y_c - 2 * y_down) / (len_c + len_down)		! From current channel surface mid-point to downstream channel surface mid-point.
			
				! Calculate flow velocity, streamflow and water volume change :
				if ( surfaceSlope .gt. 0._r8 ) then
					if ( TRunoff%wr_exchg( iu ) + ( - TRunoff%etout( iu, 1 ) - TRunoff%eroutUp( iu, 1 ) ) * Tctl%inund_dt .lt. 0._r8 ) then		! Water volume becomes negative because upward flow (from current channel to upstream channels) is too large.
						TRunoff%vr( iu, 1 ) = 0._r8
						TRunoff%erout( iu, 1 ) = 0._r8
						TRunoff%dwr( iu, 1 ) = - TRunoff%wr_exchg( iu )				
					else
						hydrR = ( TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu ) ) / ( TUnit%rwidth( iu ) + 2._r8 * TRunoff%yr_exchg( iu ) )
						TRunoff%vr( iu, 1 ) = ManningEq ( surfaceSlope, TUnit%nr( iu ), hydrR )				! Estimate flow velocity with Manning equation.
					
						TRunoff%erout( iu, 1 ) = - TRunoff%vr( iu, 1 ) * TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu )		! Q = V * A					
					
						TRunoff%dwr( iu, 1 ) = ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) + TRunoff%erout( iu, 1 ) ) * Tctl%inund_dt		! Water volume change.
					
						if ( TRunoff%erout( iu, 1 ) .lt. 0._r8 .and. TRunoff%wr_exchg( iu ) + TRunoff%dwr( iu, 1 ) .lt. 0._r8 ) then 		! Water volume becomes negative because downward flow is too large.
						
							TRunoff%erout( iu, 1 ) = - ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) + TRunoff%wr_exchg( iu ) / Tctl%inund_dt )				! All the water flows away.
							if ( TRunoff%yr_exchg( iu ) .gt. 0._r8 ) then
								TRunoff%vr( iu, 1 ) = - TRunoff%erout( iu, 1 ) / TUnit%rwidth( iu ) / TRunoff%yr_exchg( iu )
							else
								write ( *, * ) "Subroutine inund_Routing_DW ( ): Downward flow is too large & water depth is zero !"
								STOP
							end if
							TRunoff%dwr( iu, 1 ) = - TRunoff%wr_exchg( iu )
						end if
					end if
				
				elseif ( surfaceSlope .lt. 0._r8 ) then
			
					hydrR = ( TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu ) ) / ( TUnit%rwidth( iu ) + 2._r8 * TRunoff%yr_exchg( iu ) )
					TRunoff%vr( iu, 1 ) = ManningEq ( surfaceSlope, TUnit%nr( iu ), hydrR )						! Flow velocity is negative ( upward flow ).
				
					TRunoff%erout( iu, 1 ) = - TRunoff%vr( iu, 1 ) * TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu )	! Streamflow is positive ( upward flow ).				
				
					if ( TUnit%mask( iu ) .eq. 1 ) then		! This channel is NOT at basin outlet.
					
						if ( TRunoff%wr_exchg_dstrm( iu ) - TRunoff%erout( iu, 1 ) * Tctl%inund_dt .lt. 0._r8 ) then	! Upward flow (from downstream channel to current channel) is too large, so that water volume becomes negative in downstream channel.
						
							TRunoff%erout( iu, 1 ) = TRunoff%wr_exchg_dstrm( iu ) / Tctl%inund_dt
							if ( TRunoff%yr_exchg( iu ) .gt. 0._r8 ) then
								TRunoff%vr( iu, 1 ) = - TRunoff%erout( iu, 1 ) / TUnit%rwidth( iu ) / TRunoff%yr_exchg( iu )
							else
								write ( *, * ) "Subroutine inund_Routing_DW ( ): Upward flow is too large & water depth is zero !"
								STOP
							end if
						
							!TRunoff%wr_exchg_dstrm( iu ) = 0._r8			! Because all downstream-channel water flows to the current channel.
						end if
					end if
				
					TRunoff%dwr( iu, 1 ) = ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) + TRunoff%erout( iu, 1 ) ) * Tctl%inund_dt		! Water volume change.
				
				else 	! Water surface slope (surfaceSlope) is zero.
					TRunoff%vr( iu, 1 ) = 0._r8
					TRunoff%erout( iu, 1 ) = 0._r8
					TRunoff%dwr( iu, 1 ) = ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) ) * Tctl%inund_dt
				end if
			
				TRunoff%wr_rtg( iu ) = TRunoff%wr_exchg( iu ) + TRunoff%dwr( iu, 1 )						! Channel water volume after routing computation.
				TRunoff%yr_rtg( iu ) = TRunoff%wr_rtg( iu ) / TUnit%rwidth( iu ) / TUnit%rlen( iu )		! Channel water depth after routing computation.
		
			end if 		! if ( TUnit%mask( iu ) .gt. 0 )
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		
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
		integer :: k
		real( r8 ) :: hydrR						! Hydraulic radius (m).
		
		!$OMP PARALLEL FIRSTPRIVATE( iu, k, hydrR )
		!$OMP DO SCHEDULE( GUIDED )
		do iu = rtmCTL%begr, rtmCTL%endr
			if ( TUnit%mask( iu ) .gt. 0 ) then
			
				! Calculate flow velocity, streamflow and water volume change :
				hydrR = ( TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu ) ) / ( TUnit%rwidth( iu ) + 2._r8 * TRunoff%yr_exchg( iu ) )
				TRunoff%vr( iu, 1 ) = ManningEq ( TUnit%rslp( iu ), TUnit%nr( iu ), hydrR )			! Estimate flow velocity with Manning equation.
			
				TRunoff%erout( iu, 1 ) = - TRunoff%vr( iu, 1 ) * TUnit%rwidth( iu ) * TRunoff%yr_exchg( iu )			! Q = V * A			
			
				TRunoff%dwr( iu, 1 ) = ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) + TRunoff%erout( iu, 1 ) ) * Tctl%inund_dt		! Water volume change.
			
				if ( TRunoff%erout( iu, 1 ) .lt. 0._r8 .and. TRunoff%wr_exchg( iu ) + TRunoff%dwr( iu, 1 ) .lt. 0._r8 ) then 		! Water volume becomes negative because downward flow is too large.
					TRunoff%erout( iu, 1 ) = - ( - TRunoff%eroutUp( iu, 1 ) - TRunoff%etout( iu, 1 ) + TRunoff%wr_exchg( iu ) / Tctl%inund_dt )
					if ( TRunoff%yr_exchg( iu ) .gt. 0._r8 ) then
						TRunoff%vr( iu, 1 ) = - TRunoff%erout( iu, 1 ) / TUnit%rwidth( iu ) / TRunoff%yr_exchg( iu )
					else
						write ( *, * ) "Subroutine inund_Routing_KW ( ): Downward flow is too large & water depth is zero !"
						STOP
					end if
					TRunoff%dwr( iu, 1 ) = - TRunoff%wr_exchg( iu )
				end if
			
				TRunoff%wr_rtg( iu ) = TRunoff%wr_exchg( iu ) + TRunoff%dwr( iu, 1 )		! Channel water volume after routing computation.
				TRunoff%yr_rtg( iu ) = TRunoff%wr_rtg( iu ) / TUnit%rwidth( iu ) / TUnit%rlen( iu )		! Channel water depth after routing computation.
		
			end if		! if ( TUnit%mask( iu ) .gt. 0 )
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		
	end subroutine inund_Routing_KW
	
	subroutine transition2nextStep ( )
		! DESCRIPTION: Transition to the next step ( transfer values to variables of the next step ).
		implicit none
		
		TRunoff%wr( :, 1 ) = TRunoff%wr_rtg( : )			! Channel water volume.
		TRunoff%yr( :, 1 ) = TRunoff%yr_rtg( : )			! Channel water depth.
		!TRunoff%chnlQ_prev = TRunoff%chnlQ		! Channel streamflow.
		
		if ( Tctl%OPT_inund .eq. 1 ) then
			TRunoff%wf_ini = TRunoff%wf_exchg			! Floodplain water volume.
			TRunoff%hf_ini = TRunoff%hf_exchg			! Floodplain max water depth.
		end if		
	
	end subroutine transition2nextStep
	
end MODULE MOSARTinund_Core_MOD
