! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2009, Lawrence Livemore National Security Limited Liability
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! May 2015 - D. Swales - Modified for COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_ICARUS
  USE COSP_KINDS,          ONLY: wp
  USE COSP_PHYS_CONSTANTS, ONLY: amd,amw,avo,grav
  use MOD_COSP_STATS,      ONLY: hist2D
  USE MOD_COSP_CONFIG,     ONLY: R_UNDEF,numISCCPTauBins,numISCCPPresBins,isccp_histTau, &
                                 isccp_histPres
  implicit none
  
  ! Shared Parameters                   
  integer,parameter :: &
       ncolprint = 0 ! Flag for debug printing (set as parameter to increase performance)

  ! Cloud-top height determination
  integer :: &
       isccp_top_height,          & ! Top height adjustment method
       isccp_top_height_direction   ! Direction for finding atmosphere pressure level

  ! Parameters used by icarus
  real(wp),parameter :: &
       tauchk = -1._wp*log(0.9999999_wp), & ! Lower limit on optical depth
       isccp_taumin = 0.3_wp,             & ! Minimum optical depth for joint-hostogram
       pstd = 1013250._wp,                & ! Mean sea-level pressure (Pa)
       isccp_t0 = 296._wp,                & ! Mean surface temperature (K)
       output_missing_value = -1.E+30       ! Missing values

contains
  ! ##########################################################################
  ! ##########################################################################
  SUBROUTINE ICARUS(debug,debugcol,npoints,sunlit,nlev,ncol,pfull,          &
                    phalf,qv,cc,conv,dtau_s,dtau_c,th,thd,frac_out,skt,emsfc_lw,at,&
                    dem_s,dem_c,fq_isccp,totalcldarea, meanptop,meantaucld, &
                    meanalbedocld, meantb,meantbclr,boxtau,boxptop,levmatch)
        
    ! INPUTS 
    INTEGER,intent(in) ::      & !
         npoints,              & ! Number of model points in the horizontal
         nlev,                 & ! Number of model levels in column
         ncol,                 & ! Number of subcolumns
         debug,                & ! Debug flag
         debugcol                ! Debug column flag
    INTEGER,intent(in),dimension(npoints) :: & !
         sunlit                  ! 1 for day points, 0 for night time 
    REAL(WP),intent(in) ::     & !
         emsfc_lw                ! 10.5 micron emissivity of surface (fraction)  
    REAL(WP),intent(in),dimension(npoints) :: & !
         skt                     ! Skin Temperature (K)
    REAL(WP),intent(in),dimension(npoints,ncol,nlev) :: & !
         frac_out                ! Boxes gridbox divided up into subcolumns
    REAL(WP),intent(in),dimension(npoints,nlev) :: & !
         pfull,                & ! Pressure of full model levels (Pascals)
         qv,                   & ! Water vapor specific humidity (kg vapor/ kg air)
         cc,                   & ! Cloud cover in each model level (fraction)
         conv,                 & ! Convective cloud cover in each model
         at,                   & ! Temperature in each model level (K) 
         dem_c,                & ! Emissivity for convective clouds
         dem_s,                & ! Emissivity for stratiform clouds
         dtau_c,               & ! Optical depth for convective clouds
         dtau_s                  ! Optical depth for stratiform clouds 
    REAL(WP),intent(in),dimension(npoints,nlev+1) :: & !
         phalf                   ! Pressure of half model levels (Pascals)!
    integer,intent(in) :: th,thd

    ! OUTPUTS 
    REAL(WP),intent(out),dimension(npoints,7,7) :: & 
         fq_isccp                ! The fraction of the model grid box covered by clouds
    REAL(WP),intent(out),dimension(npoints) :: & 
         totalcldarea,         & ! The fraction of model grid box columns with cloud present
         meanptop,             & ! Mean cloud top pressure (mb) - linear averaging
         meantaucld,           & ! Mean optical thickness
         meanalbedocld,        & ! Mean cloud albedo  
         meantb,               & ! Mean all-sky 10.5 micron brightness temperature
         meantbclr               ! Mean clear-sky 10.5 micron brightness temperature
    REAL(WP),intent(out),dimension(npoints,ncol) :: & 
         boxtau,               & ! Optical thickness in each column
         boxptop                 ! Cloud top pressure (mb) in each column
    INTEGER,intent(out),dimension(npoints,ncol) :: &
         levmatch                ! Used for icarus unit testing only


    ! INTERNAL VARIABLES
    CHARACTER(len=10)                     :: ftn09
    REAL(WP),dimension(npoints,ncol)      :: boxttop
    REAL(WP),dimension(npoints,ncol,nlev) :: dtau,demIN
    INTEGER                               :: j,ilev,ibox
    INTEGER,dimension(nlev,ncol   )       :: acc

    ! PARAMETERS
    character ,parameter, dimension(6) :: cchar=(/' ','-','1','+','I','+'/)
    character(len=1),parameter,dimension(6) :: cchar_realtops=(/ ' ',' ','1','1','I','I'/)
    ! ##########################################################################
    
    call cosp_simulator_optics(npoints,ncol,nlev,frac_out,dem_c,dem_s,demIN)
    call cosp_simulator_optics(npoints,ncol,nlev,frac_out,dtau_c,dtau_s,dtau)

    call ICARUS_SUBCOLUMN(npoints,ncol,nlev,sunlit,dtau,demIN,skt,emsfc_lw,qv,at,        &
                      pfull,phalf,frac_out,levmatch,boxtau,boxptop,boxttop,meantbclr)

    call ICARUS_COLUMN(npoints,ncol,boxtau,boxptop/100._wp,sunlit,boxttop,&
                       fq_isccp,meanalbedocld,meanptop,meantaucld,totalcldarea,meantb)

    ! ##########################################################################
    ! OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
    ! ##########################################################################
    
    if (debugcol.ne.0) then
       do j=1,npoints,debugcol
          
          ! Produce character output
          do ilev=1,nlev
             acc(ilev,1:ncol)=frac_out(j,1:ncol,ilev)*2
             where(levmatch(j,1:ncol) .eq. ilev) acc(ilev,1:ncol)=acc(ilev,1:ncol)+1
          enddo
          
          write(ftn09,11) j
11        format('ftn09.',i4.4)
          open(9, FILE=ftn09, FORM='FORMATTED')
          
          write(9,'(a1)') ' '
          write(9,'(10i5)') (ilev,ilev=5,nlev,5)
          write(9,'(a1)') ' '
          
          do ibox=1,ncol
             write(9,'(40(a1),1x,40(a1))') &
                  (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev),&
                  (cchar(acc(ilev,ibox)+1),ilev=1,nlev) 
          end do
          close(9)

       enddo       
    end if
    
    return
  end SUBROUTINE ICARUS
  
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  SUBROUTINE ICARUS_SUBCOLUMN(npoints,ncol,nlev,sunlit,dtau,demiN,skt,emsfc_lw,qv,at,        &
                          pfull,phalf,frac_out,levmatch,boxtau,boxptop,boxttop,meantbclr)
    ! Inputs
    INTEGER, intent(in) ::   &
         ncol,               & ! Number of subcolumns
         npoints,            & ! Number of horizontal gridpoints
         nlev                  ! Number of vertical levels
    INTEGER, intent(in), dimension(npoints) :: &
         sunlit                ! 1=day 0=night
    REAL(WP),intent(in) :: &
         emsfc_lw              ! 10.5 micron emissivity of surface (fraction) 
    REAL(WP),intent(in), dimension(npoints) ::  &
         skt                   ! Skin temperature
    REAL(WP),intent(in), dimension(npoints,nlev) ::  &
         at,                 & ! Temperature 
         pfull,              & ! Presure
         qv                    ! Specific humidity
    REAL(WP),intent(in), dimension(npoints,ncol,nlev) :: &
         frac_out,           & ! Subcolumn cloud cover
         dtau,               & ! Subcolumn optical thickness
         demIN                 ! Subcolumn emissivity
    REAL(WP),intent(in), dimension(npoints,nlev+1) :: &
         phalf                 ! Pressure at model half levels

    ! Outputs
    REAL(WP),intent(inout),dimension(npoints) :: &
         meantbclr             ! Mean clear-sky 10.5 micron brightness temperature
    REAL(WP),intent(inout),dimension(npoints,ncol) :: &
         boxtau,             & ! Optical thickness in each column
         boxptop,            & ! Cloud top pressure (mb) in each column
         boxttop               ! Cloud top temperature in each column
    INTEGER, intent(inout),dimension(npoints,ncol)      :: levmatch

    ! Local Variables
    INTEGER :: &
       j,ibox,ilev,k1,k2,icycle
    INTEGER,dimension(npoints) :: &
       nmatch,itrop
    INTEGER,dimension(npoints,nlev-1) :: &
       match
    REAL(WP) :: &
       logp,logp1,logp2,atd
    REAL(WP),dimension(npoints) :: &
       bb,attropmin,attrop,ptrop,atmax,btcmin,transmax,tauir,taumin,fluxtopinit,press,   &
       dpress,atmden,rvh20,rhoave,rh20s,rfrgn,tmpexp,tauwv,wk,trans_layers_above_clrsky, &
       fluxtop_clrsky
    REAL(WP),dimension(npoints,nlev) :: &
       dem_wv
    REAL(WP),dimension(npoints,ncol) :: &
       trans_layers_above,dem,tb,emcld,fluxtop,tau,ptop

    ! ####################################################################################
    ! Compute cloud optical depth for each column by summing up subcolumns
    tau(1:npoints,1:ncol) = 0._wp
    tau(1:npoints,1:ncol) = sum(dtau,dim=3)

    ! Set tropopause values
    if (isccp_top_height .eq. 1 .or. isccp_top_height .eq. 3) then 
       ptrop(1:npoints)     = 5000._wp
       attropmin(1:npoints) = 400._wp
       atmax(1:npoints)     = 0._wp
       attrop(1:npoints)    = 120._wp
       itrop(1:npoints)     = 1

       do ilev=1,nlev
          where(pfull(1:npoints,ilev) .lt. 40000. .and. &
                pfull(1:npoints,ilev) .gt.  5000. .and. &
                at(1:npoints,ilev)    .lt. attropmin(1:npoints))
             ptrop(1:npoints)     = pfull(1:npoints,ilev)
             attropmin(1:npoints) = at(1:npoints,ilev)
             attrop(1:npoints)    = attropmin(1:npoints)
             itrop     = ilev
          endwhere
       enddo

       do ilev=1,nlev
          atmax(1:npoints) = merge(at(1:npoints,ilev),atmax(1:npoints),&
               at(1:npoints,ilev) .gt. atmax(1:npoints) .and. ilev  .ge. itrop(1:npoints))
       enddo
    end if
  
    if (isccp_top_height .eq. 1 .or. isccp_top_height .eq. 3) then
       ! ############################################################################
       !                        Clear-sky radiance calculation
       !       
       ! Compute water vapor continuum emissivity this treatment follows Schwarkzopf 
       ! and Ramasamy JGR 1999,vol 104, pages 9467-9499. The emissivity is calculated
       ! at a wavenumber of 955 cm-1, or 10.47 microns 
       ! ############################################################################
       do ilev=1,nlev
          press(1:npoints)  = pfull(1:npoints,ilev)*10._wp
          dpress(1:npoints) = (phalf(1:npoints,ilev+1)-phalf(1:npoints,ilev))*10
          atmden(1:npoints) = dpress(1:npoints)/(grav*100._wp)
          rvh20(1:npoints)  = qv(1:npoints,ilev)*amd/amw
          wk(1:npoints)     = rvh20(1:npoints)*avo*atmden/amd
          rhoave(1:npoints) = (press(1:npoints)/pstd)*(isccp_t0/at(1:npoints,ilev))
          rh20s(1:npoints)  = rvh20(1:npoints)*rhoave(1:npoints)
          rfrgn(1:npoints)  = rhoave(1:npoints)-rh20s(1:npoints)
          tmpexp(1:npoints) = exp(-0.02_wp*(at(1:npoints,ilev)-isccp_t0))
          tauwv(1:npoints)  = wk(1:npoints)*1.e-20*((0.0224697_wp*rh20s(1:npoints)*      &
                              tmpexp(1:npoints))+(3.41817e-7*rfrgn(1:npoints)))*0.98_wp
          dem_wv(1:npoints,ilev) = 1._wp - exp( -1._wp * tauwv(1:npoints))
       enddo

       fluxtop_clrsky(1:npoints)            = 0._wp
       trans_layers_above_clrsky(1:npoints) = 1._wp
       do ilev=1,nlev
          ! Black body emission at temperature of the layer
          bb(1:npoints) = 1._wp / ( exp(1307.27_wp/at(1:npoints,ilev)) - 1._wp )
          
          ! Increase TOA flux by flux emitted from layer times total transmittance in layers above
          fluxtop_clrsky(1:npoints) = fluxtop_clrsky(1:npoints) + &
               dem_wv(1:npoints,ilev)*bb(1:npoints)*trans_layers_above_clrsky(1:npoints)
          
          ! Update trans_layers_above with transmissivity from this layer for next time around loop
          trans_layers_above_clrsky(1:npoints) = trans_layers_above_clrsky(1:npoints)*&
              (1.-dem_wv(1:npoints,ilev))                
       enddo

       ! Add in surface emission
       bb(1:npoints) = 1._wp/( exp(1307.27_wp/skt(1:npoints)) - 1._wp )
       fluxtop_clrsky(1:npoints) = fluxtop_clrsky(1:npoints) + &
           emsfc_lw * bb(1:npoints)*trans_layers_above_clrsky(1:npoints)

       ! Clear Sky brightness temperature
       meantbclr(1:npoints) = 1307.27_wp/(log(1._wp+(1._wp/fluxtop_clrsky(1:npoints))))
       
       ! #################################################################################
       !                        All-sky radiance calculation
       ! #################################################################################
       
       fluxtop(1:npoints,1:ncol)            = 0._wp
       trans_layers_above(1:npoints,1:ncol) = 1._wp
       do ilev=1,nlev
          ! Black body emission at temperature of the layer
          bb=1._wp/(exp(1307.27_wp/at(1:npoints,ilev)) - 1._wp)
          
          do ibox=1,ncol
             ! Emissivity
             dem(1:npoints,ibox) = merge(dem_wv(1:npoints,ilev), &
                                         1._wp-(1._wp-demIN(1:npoints,ibox,ilev))*(1._wp-dem_wv(1:npoints,ilev)), &
                                         demIN(1:npoints,ibox,ilev) .eq. 0)

             ! Increase TOA flux emitted from layer
             fluxtop(1:npoints,ibox) = fluxtop(1:npoints,ibox) + dem(1:npoints,ibox)*bb*trans_layers_above(1:npoints,ibox) 
             
             ! Update trans_layer by emitted layer from above
             trans_layers_above(1:npoints,ibox) = trans_layers_above(1:npoints,ibox)*(1._wp-dem(1:npoints,ibox))
          enddo
       enddo

       ! Add in surface emission
       bb(1:npoints)=1._wp/( exp(1307.27_wp/skt(1:npoints)) - 1._wp )
       do ibox=1,ncol
          fluxtop(1:npoints,ibox) = fluxtop(1:npoints,ibox) + emsfc_lw*bb(1:npoints)*trans_layers_above(1:npoints,ibox) 
       end do

       ! All Sky brightness temperature
       boxttop(1:npoints,1:ncol) = 1307.27_wp/(log(1._wp+(1._wp/fluxtop(1:npoints,1:ncol))))

       ! #################################################################################  
       !                            Cloud-Top Temperature
       !
       ! Now that you have the top of atmosphere radiance, account for ISCCP 
       ! procedures to determine cloud top temperature account for partially
       ! transmitting cloud recompute flux ISCCP would see assuming a single layer
       ! cloud. *NOTE* choice here of 2.13, as it is primarily ice clouds which have 
       ! partial emissivity and need the adjustment performed in this section. If it
       ! turns out that the cloud brightness temperature is greater than 260K, then 
       ! the liquid cloud conversion factor of 2.56 is used. *NOTE* that this is 
       ! discussed on pages 85-87 of the ISCCP D level documentation 
       ! (Rossow et al. 1996)
       ! #################################################################################

       ! Compute minimum brightness temperature and optical depth
       btcmin(1:npoints) = 1._wp /  ( exp(1307.27_wp/(attrop(1:npoints)-5._wp)) - 1._wp ) 

       do ibox=1,ncol
          transmax(1:npoints) = (fluxtop(1:npoints,ibox)-btcmin) /(fluxtop_clrsky(1:npoints)-btcmin(1:npoints))
          tauir(1:npoints)    = tau(1:npoints,ibox)/2.13_wp
          taumin(1:npoints)   = -log(max(min(transmax(1:npoints),0.9999999_wp),0.001_wp))
          if (isccp_top_height .eq. 1) then
             do j=1,npoints  
                if (transmax(j) .gt. 0.001 .and.  transmax(j) .le. 0.9999999) then
                   fluxtopinit(j) = fluxtop(j,ibox)
                   tauir(j) = tau(j,ibox)/2.13_wp
                endif
             enddo
             do icycle=1,2
                do j=1,npoints  
                   if (tau(j,ibox) .gt. (tauchk)) then 
                      if (transmax(j) .gt. 0.001 .and.  transmax(j) .le. 0.9999999) then
                         emcld(j,ibox) = 1._wp - exp(-1._wp * tauir(j)  )
                         fluxtop(j,ibox) = fluxtopinit(j) - ((1.-emcld(j,ibox))*fluxtop_clrsky(j))
                         fluxtop(j,ibox)=max(1.E-06_wp,(fluxtop(j,ibox)/emcld(j,ibox)))
                         tb(j,ibox)= 1307.27_wp / (log(1._wp + (1._wp/fluxtop(j,ibox))))
                         if (tb(j,ibox) .gt. 260.) then
                            tauir(j) = tau(j,ibox) / 2.56_wp
                         end if
                      end if
                   end if
                enddo
            enddo
          endif

          ! Cloud-top temperature
          where(tau(1:npoints,ibox) .gt. tauchk)
             tb(1:npoints,ibox)= 1307.27_wp/ (log(1. + (1._wp/fluxtop(1:npoints,ibox))))
             where (isccp_top_height .eq. 1 .and. tauir(1:npoints) .lt. taumin(1:npoints))
                tb(1:npoints,ibox) = attrop(1:npoints) - 5._wp 
                tau(1:npoints,ibox) = 2.13_wp*taumin(1:npoints)
             endwhere
          endwhere
          
          ! Clear-sky brightness temperature
          where(tau(1:npoints,ibox) .le. tauchk) 
             tb(1:npoints,ibox) = meantbclr(1:npoints)
          endwhere
       enddo
    else
       meantbclr(1:npoints) = output_missing_value
    end if

    ! ####################################################################################
    !                           Cloud-Top Pressure
    !
    ! The 2 methods differ according to whether or not you use the physical cloud
    ! top pressure (isccp_top_height = 2) or the radiatively determined cloud top
    ! pressure (isccp_top_height = 1 or 3)
    ! ####################################################################################
    do ibox=1,ncol
       !segregate according to optical thickness
       if (isccp_top_height .eq. 1 .or. isccp_top_height .eq. 3) then  
          
          ! Find level whose temperature most closely matches brightness temperature
          nmatch(1:npoints)=0
          do k1=1,nlev-1
             ilev = merge(nlev-k1,k1,isccp_top_height_direction .eq. 2)        
             do j=1,npoints 
                if (ilev           .ge. itrop(j)     .and. &
                     ((at(j,ilev)  .ge. tb(j,ibox)   .and. &  
                      at(j,ilev+1) .le. tb(j,ibox))  .or.  &
                      (at(j,ilev)  .le. tb(j,ibox)   .and. &
                      at(j,ilev+1) .ge. tb(j,ibox)))) then 
                   nmatch(j)=nmatch(j)+1
                   match(j,nmatch(j))=ilev
                endif
             enddo
          enddo

          do j=1,npoints 
             if (nmatch(j) .ge. 1) then
                k1 = match(j,nmatch(j))
                k2 = k1 + 1
                logp1 = log(pfull(j,k1))
                logp2 = log(pfull(j,k2))
                atd = max(tauchk,abs(at(j,k2) - at(j,k1)))
                logp=logp1+(logp2-logp1)*abs(tb(j,ibox)-at(j,k1))/atd
                ptop(j,ibox) = exp(logp)
                levmatch(j,ibox) = merge(k1,k2,abs(pfull(j,k1)-ptop(j,ibox)) .lt. abs(pfull(j,k2)-ptop(j,ibox)))
             else
                if (tb(j,ibox) .le. attrop(j)) then
                   ptop(j,ibox)=ptrop(j)
                   levmatch(j,ibox)=itrop(j)
                end if
                if (tb(j,ibox) .ge. atmax(j)) then
                   ptop(j,ibox)=pfull(j,nlev)
                   levmatch(j,ibox)=nlev
                end if
             end if
          enddo
       else
          ptop(1:npoints,ibox)=0.
          do ilev=1,nlev
             where((ptop(1:npoints,ibox) .eq. 0. ) .and.(frac_out(1:npoints,ibox,ilev) .ne. 0))
                ptop(1:npoints,ibox)=phalf(1:npoints,ilev)
                levmatch(1:npoints,ibox)=ilev
             endwhere
          end do
       end if
       where(tau(1:npoints,ibox) .le. tauchk)
          ptop(1:npoints,ibox)=0._wp
          levmatch(1:npoints,ibox)=0._wp
       endwhere
    enddo

    ! ####################################################################################
    !                Compute subcolumn pressure and optical depth
    ! ####################################################################################
    boxtau(1:npoints,1:ncol)  = output_missing_value
    boxptop(1:npoints,1:ncol) = output_missing_value
    do ibox=1,ncol
       do j=1,npoints 
          if (tau(j,ibox) .gt. (tauchk) .and. ptop(j,ibox) .gt. 0.) then
             if (sunlit(j).eq.1 .or. isccp_top_height .eq. 3) then
                boxtau(j,ibox) = tau(j,ibox)
                boxptop(j,ibox) = ptop(j,ibox)!/100._wp
             endif
          endif
       enddo
    enddo

  end SUBROUTINE ICARUS_SUBCOLUMN

  ! ######################################################################################
  ! SUBROUTINE icarus_column
  ! ######################################################################################
  SUBROUTINE ICARUS_column(npoints,ncol,boxtau,boxptop,sunlit,boxttop,fq_isccp,     &
                           meanalbedocld,meanptop,meantaucld,totalcldarea,meantb)
    ! Inputs
    INTEGER, intent(in) :: &
         ncol,    & ! Number of subcolumns
         npoints    ! Number of horizontal gridpoints
    INTEGER, intent(in),dimension(npoints) :: &
         sunlit     ! day=1 night=0 
    REAL(WP),intent(in),dimension(npoints,ncol) ::  &
         boxttop,  & ! Subcolumn top temperature
         boxptop,  & ! Subcolumn cloud top pressure
         boxtau      ! Subcolumn optical depth

    ! Outputs
    REAL(WP),intent(inout),dimension(npoints) :: &
         meanalbedocld, & ! Gridmean cloud albedo
         meanptop,      & ! Gridmean cloud top pressure (mb) - linear averaging
         meantaucld,    & ! Gridmean optical thickness
         totalcldarea,  & ! The fraction of model grid box columns with cloud present
         meantb           ! Gridmean all-sky 10.5 micron brightness temperature 
    REAL(WP),intent(inout),dimension(npoints,7,7) :: &
         fq_isccp         ! The fraction of the model grid box covered by clouds

    ! Local Variables
    INTEGER :: j,ilev,ilev2
    REAL(WP),dimension(npoints,ncol) :: albedocld
    LOGICAL, dimension(npoints,ncol) :: box_cloudy

    ! Variables for new joint-histogram implementation
    logical,dimension(ncol) :: box_cloudy2

    ! ####################################################################################
    !                           Brightness Temperature
    ! ####################################################################################
    if (isccp_top_height .eq. 1 .or. isccp_top_height .eq. 3) then
       meantb(1:npoints)=sum(boxttop,2)/ncol
    else
       meantb(1:npoints) = output_missing_value
    endif

    ! ####################################################################################
    !                 Determines ISCCP cloud type frequencies
    !
    ! Now that boxptop and boxtau have been determined, determine amount of each of the 
    ! 49 ISCCP cloud types. Also compute grid box mean cloud top pressure and 
    ! optical thickness.  The mean cloud top pressure and optical thickness are 
    ! averages over the cloudy area only. The mean cloud top pressure is a linear
    ! average of the cloud top pressures. The mean cloud optical thickness is 
    ! computed by converting optical thickness to an albedo, averaging in albedo 
    ! units, then converting the average albedo back to a mean optical thickness.  
    ! ####################################################################################

    ! Initialize
    albedocld(1:npoints,1:ncol)  = 0._wp
    box_cloudy(1:npoints,1:ncol) = .false.
    
    ! Reset frequencies
    !fq_isccp = spread(spread(merge(0._wp,output_missing_value,sunlit .eq. 1 .or. isccp_top_height .eq. 3),2,7),2,7)
    do ilev=1,7
       do ilev2=1,7
          do j=1,npoints ! 
             if (sunlit(j).eq.1 .or. isccp_top_height .eq. 3) then 
                fq_isccp(j,ilev,ilev2)= 0.
	     else 
                fq_isccp(j,ilev,ilev2)= output_missing_value
             end if
          enddo
       enddo
    enddo

    
    ! Reset variables need for averaging cloud properties
    where(sunlit .eq. 1 .or. isccp_top_height .eq. 3)
       totalcldarea(1:npoints)  = 0._wp
       meanalbedocld(1:npoints) = 0._wp
       meanptop(1:npoints)      = 0._wp
       meantaucld(1:npoints)    = 0._wp
    elsewhere
       totalcldarea(1:npoints)  = output_missing_value
       meanalbedocld(1:npoints) = output_missing_value
       meanptop(1:npoints)      = output_missing_value
       meantaucld(1:npoints)    = output_missing_value
    endwhere
    
    ! Compute column quantities and joint-histogram
    do j=1,npoints 
       ! Subcolumns that are cloudy(true) and not(false)
       box_cloudy2(1:ncol) = merge(.true.,.false.,boxtau(j,1:ncol) .gt. tauchk .and. boxptop(j,1:ncol) .gt. 0.)

       ! Compute joint histogram and column quantities for points that are sunlit and cloudy
       if (sunlit(j) .eq.1 .or. isccp_top_height .eq. 3) then 
          ! Joint-histogram
          call hist2D(boxtau(j,1:ncol),boxptop(j,1:ncol),ncol,isccp_histTau,numISCCPTauBins, &
               isccp_histPres,numISCCPPresBins,fq_isccp(j,1:numISCCPTauBins,1:numISCCPPresBins))
          fq_isccp(j,1:numISCCPTauBins,1:numISCCPPresBins) = &
               fq_isccp(j,1:numISCCPTauBins,1:numISCCPPresBins)/ncol
          
          ! Column cloud area
          totalcldarea(j) = real(count(box_cloudy2(1:ncol) .and. boxtau(j,1:ncol) .gt. isccp_taumin))/ncol
             
          ! Subcolumn cloud albedo
          !albedocld(j,1:ncol) = merge((boxtau(j,1:ncol)**0.895_wp)/((boxtau(j,1:ncol)**0.895_wp)+6.82_wp),&
          !     0._wp,box_cloudy2(1:ncol) .and. boxtau(j,1:ncol) .gt. isccp_taumin)
          where(box_cloudy2(1:ncol) .and. boxtau(j,1:ncol) .gt. isccp_taumin)
             albedocld(j,1:ncol) = (boxtau(j,1:ncol)**0.895_wp)/((boxtau(j,1:ncol)**0.895_wp)+6.82_wp)
          elsewhere
             albedocld(j,1:ncol) = 0._wp
          endwhere
          
          ! Column cloud albedo
          meanalbedocld(j) = sum(albedocld(j,1:ncol))/ncol
          
          ! Column cloud top pressure
          meanptop(j) = sum(boxptop(j,1:ncol),box_cloudy2(1:ncol) .and. boxtau(j,1:ncol) .gt. isccp_taumin)/ncol
       endif
    enddo
    
    ! Compute mean cloud properties. Set to mssing value in the event that totalcldarea=0
    where(totalcldarea(1:npoints) .gt. 0)
       meanptop(1:npoints)      = 100._wp*meanptop(1:npoints)/totalcldarea(1:npoints)
       meanalbedocld(1:npoints) = meanalbedocld(1:npoints)/totalcldarea(1:npoints)
       meantaucld(1:npoints)    = (6.82_wp/((1._wp/meanalbedocld(1:npoints))-1.))**(1._wp/0.895_wp)
    elsewhere
       meanptop(1:nPoints)      = output_missing_value
       meanalbedocld(1:nPoints) = output_missing_value
       meantaucld(1:nPoints)    = output_missing_value
    endwhere
    !meanptop(1:npoints)      = merge(100._wp*meanptop(1:npoints)/totalcldarea(1:npoints),&
    !                                 output_missing_value,totalcldarea(1:npoints) .gt. 0)
    !meanalbedocld(1:npoints) = merge(meanalbedocld(1:npoints)/totalcldarea(1:npoints), &
    !                                 output_missing_value,totalcldarea(1:npoints) .gt. 0)
    !meantaucld(1:npoints)    = merge((6.82_wp/((1._wp/meanalbedocld(1:npoints))-1.))**(1._wp/0.895_wp), &
    !                                 output_missing_value,totalcldarea(1:npoints) .gt. 0)

    ! Represent in percent
    where(totalcldarea .ne. output_missing_value) totalcldarea = totalcldarea*100._wp
    where(fq_isccp     .ne. output_missing_value) fq_isccp     = fq_isccp*100._wp
    
    
  end SUBROUTINE ICARUS_column
  
  subroutine cosp_simulator_optics(dim1,dim2,dim3,flag,varIN1,varIN2,varOUT)
    ! INPUTS
    integer,intent(in) :: &
         dim1,   & ! Dimension 1 extent (Horizontal)
         dim2,   & ! Dimension 2 extent (Subcolumn)
         dim3      ! Dimension 3 extent (Vertical)
    real(wp),intent(in),dimension(dim1,dim2,dim3) :: &
         flag      ! Logical to determine the of merge var1IN and var2IN
    real(wp),intent(in),dimension(dim1,     dim3) :: &
         varIN1, & ! Input field 1
         varIN2    ! Input field 2
    ! OUTPUTS
    real(wp),intent(out),dimension(dim1,dim2,dim3) :: &
         varOUT    ! Merged output field
    ! LOCAL VARIABLES
    integer :: j
    
    varOUT(1:dim1,1:dim2,1:dim3) = 0._wp
    do j=1,dim2
       where(flag(:,j,:) .eq. 1)
          varOUT(:,j,:) = varIN2
       endwhere
       where(flag(:,j,:) .eq. 2)
          varOUT(:,j,:) = varIN1
       endwhere
    enddo
  end subroutine cosp_simulator_optics
end module MOD_ICARUS

