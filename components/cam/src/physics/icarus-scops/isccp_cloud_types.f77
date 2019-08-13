      SUBROUTINE ISCCP_CLOUD_TYPES(
     &     debug,
     &     debugcol,
     &     npoints,
     &     sunlit,
     &     nlev,
     &     ncol,
     &     seed,
     &     pfull,
     &     phalf,
     &     qv,
     &     cc,
     &     conv,
     &     dtau_s,
     &     dtau_c,
     &     top_height,
     &     top_height_direction,
     &     overlap,
     &     frac_out,
     &     skt,
     &     emsfc_lw,
     &     at,
     &     dem_s,
     &     dem_c,
     &     fq_isccp,
     &     totalcldarea,
     &     meanptop,
     &     meantaucld,
     &     meanalbedocld,
     &     meantb,
     &     meantbclr,
     &     boxtau,
     &     boxptop
     &)

!$Id: isccp_cloud_types.f,v 4.0 2009/03/06 11:05:11 hadmw Exp $

! *****************************COPYRIGHT****************************
! (c) British Crown Copyright 2009, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the
! following conditions are met:
! 
!     * Redistributions of source code must retain the above 
!       copyright  notice, this list of conditions and the following 
!       disclaimer.
!     * Redistributions in binary form must reproduce the above 
!       copyright notice, this list of conditions and the following 
!       disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its 
!       contributors may be used to endorse or promote products
!       derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
! 
! *****************************COPYRIGHT*******************************
! *****************************COPYRIGHT*******************************
! *****************************COPYRIGHT*******************************

      implicit none

!     NOTE:   the maximum number of levels and columns is set by
!             the following parameter statement

      INTEGER ncolprint
      
!     -----
!     Input 
!     -----

      INTEGER npoints       !  number of model points in the horizontal
      INTEGER nlev          !  number of model levels in column
      INTEGER ncol          !  number of subcolumns

      INTEGER sunlit(npoints) !  1 for day points, 0 for night time

      INTEGER seed(npoints)
      !  seed values for marsaglia  random number generator
      !  It is recommended that the seed is set
      !  to a different value for each model
      !  gridbox it is called on, as it is
      !  possible that the choice of the same
      !  seed value every time may introduce some
      !  statistical bias in the results, particularly
      !  for low values of NCOL.

      REAL pfull(npoints,nlev)
                       !  pressure of full model levels (Pascals)
                  !  pfull(npoints,1) is top level of model
                  !  pfull(npoints,nlev) is bot of model

      REAL phalf(npoints,nlev+1)
                  !  pressure of half model levels (Pascals)
                  !  phalf(npoints,1) is top of model
                  !  phalf(npoints,nlev+1) is the surface pressure

      REAL qv(npoints,nlev)
                  !  water vapor specific humidity (kg vapor/ kg air)
                  !         on full model levels

      REAL cc(npoints,nlev)   
                  !  input cloud cover in each model level (fraction) 
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by clouds

      REAL conv(npoints,nlev) 
                  !  input convective cloud cover in each model
                  !   level (fraction) 
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by convective clouds

      REAL dtau_s(npoints,nlev) 
                  !  mean 0.67 micron optical depth of stratiform
                !  clouds in each model level
                  !  NOTE:  this the cloud optical depth of only the
                  !  cloudy part of the grid box, it is not weighted
                  !  with the 0 cloud optical depth of the clear
                  !         part of the grid box

      REAL dtau_c(npoints,nlev) 
                  !  mean 0.67 micron optical depth of convective
                !  clouds in each
                  !  model level.  Same note applies as in dtau_s.

      INTEGER overlap                   !  overlap type
                              !  1=max
                              !  2=rand
                              !  3=max/rand

      INTEGER top_height                !  1 = adjust top height using both a computed
                                        !  infrared brightness temperature and the visible
                              !  optical depth to adjust cloud top pressure. Note
                              !  that this calculation is most appropriate to compare
                              !  to ISCCP data during sunlit hours.
                                        !  2 = do not adjust top height, that is cloud top
                                        !  pressure is the actual cloud top pressure
                                        !  in the model
                              !  3 = adjust top height using only the computed
                              !  infrared brightness temperature. Note that this
                              !  calculation is most appropriate to compare to ISCCP
                              !  IR only algortihm (i.e. you can compare to nighttime
                              !  ISCCP data with this option)

      INTEGER top_height_direction ! direction for finding atmosphere pressure level
                                 ! with interpolated temperature equal to the radiance
				 ! determined cloud-top temperature
				 !
				 ! 1 = find the *lowest* altitude (highest pressure) level
				 ! with interpolated temperature equal to the radiance
				 ! determined cloud-top temperature
				 !
				 ! 2 = find the *highest* altitude (lowest pressure) level
				 ! with interpolated temperature equal to the radiance 
				 ! determined cloud-top temperature
				 ! 
				 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
				 !
				 ! 1 = old setting: matches all versions of 
				 ! ISCCP simulator with versions numbers 3.5.1 and lower
				 !
				 ! 2 = default setting: for version numbers 4.0 and higher  
!
!     The following input variables are used only if top_height = 1 or top_height = 3
!
      REAL skt(npoints)                 !  skin Temperature (K)
      REAL emsfc_lw                     !  10.5 micron emissivity of surface (fraction)                                            
      REAL at(npoints,nlev)                   !  temperature in each model level (K)
      REAL dem_s(npoints,nlev)                !  10.5 micron longwave emissivity of stratiform
                              !  clouds in each
                                        !  model level.  Same note applies as in dtau_s.
      REAL dem_c(npoints,nlev)                  !  10.5 micron longwave emissivity of convective
                              !  clouds in each
                                        !  model level.  Same note applies as in dtau_s.

      REAL frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                              ! Equivalent of BOX in original version, but
                              ! indexed by column then row, rather than
                              ! by row then column



!     ------
!     Output
!     ------

      REAL fq_isccp(npoints,7,7)        !  the fraction of the model grid box covered by
                                        !  each of the 49 ISCCP D level cloud types

      REAL totalcldarea(npoints)        !  the fraction of model grid box columns
                                        !  with cloud somewhere in them.  NOTE: This diagnostic
					! does not count model clouds with tau < isccp_taumin
                              ! Thus this diagnostic does not equal the sum over all entries of fq_isccp.
			      ! However, this diagnostic does equal the sum over entries of fq_isccp with
			      ! itau = 2:7 (omitting itau = 1)
      
      
      ! The following three means are averages only over the cloudy areas with tau > isccp_taumin.  
      ! If no clouds with tau > isccp_taumin are in grid box all three quantities should equal zero.      
                              
      REAL meanptop(npoints)            !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
                              
      REAL meantaucld(npoints)          !  mean optical thickness 
                                        !  linear averaging in albedo performed.
      
      real meanalbedocld(npoints)        ! mean cloud albedo
                                        ! linear averaging in albedo performed
					
      real meantb(npoints)              ! mean all-sky 10.5 micron brightness temperature
      
      real meantbclr(npoints)           ! mean clear-sky 10.5 micron brightness temperature
      
      REAL boxtau(npoints,ncol)         !  optical thickness in each column
      
      REAL boxptop(npoints,ncol)        !  cloud top pressure (mb) in each column
                              
                                                                                          
!
!     ------
!     Working variables added when program updated to mimic Mark Webb's PV-Wave code
!     ------

      REAL dem(npoints,ncol),bb(npoints)     !  working variables for 10.5 micron longwave 
                              !  emissivity in part of
                              !  gridbox under consideration

      REAL ptrop(npoints)
      REAL attrop(npoints)
      REAL attropmin (npoints)
      REAL atmax(npoints)
      REAL atmin(npoints)
      REAL btcmin(npoints)
      REAL transmax(npoints)

      INTEGER i,j,ilev,ibox,itrop(npoints)
      INTEGER ipres(npoints)
      INTEGER itau(npoints),ilev2
      INTEGER acc(nlev,ncol)
      INTEGER match(npoints,nlev-1)
      INTEGER nmatch(npoints)
      INTEGER levmatch(npoints,ncol)
      
      !variables needed for water vapor continuum absorption
      real fluxtop_clrsky(npoints),trans_layers_above_clrsky(npoints)
      real taumin(npoints)
      real dem_wv(npoints,nlev), wtmair, wtmh20, Navo, grav, pstd, t0
      real press(npoints), dpress(npoints), atmden(npoints)
      real rvh20(npoints), wk(npoints), rhoave(npoints)
      real rh20s(npoints), rfrgn(npoints)
      real tmpexp(npoints),tauwv(npoints)
      
      character*1 cchar(6),cchar_realtops(6)
      integer icycle
      REAL tau(npoints,ncol)
      LOGICAL box_cloudy(npoints,ncol)
      REAL tb(npoints,ncol)
      REAL ptop(npoints,ncol)
      REAL emcld(npoints,ncol)
      REAL fluxtop(npoints,ncol)
      REAL trans_layers_above(npoints,ncol)
      real isccp_taumin,fluxtopinit(npoints),tauir(npoints)
      REAL albedocld(npoints,ncol)
      real boxarea
      integer debug       ! set to non-zero value to print out inputs
                    ! with step debug
      integer debugcol    ! set to non-zero value to print out column
                    ! decomposition with step debugcol
      integer rangevec(npoints),rangeerror

      integer index1(npoints),num1,jj,k1,k2
      real rec2p13,tauchk,logp,logp1,logp2,atd

      character*10 ftn09
      
      DATA isccp_taumin / 0.3 /
      DATA cchar / ' ','-','1','+','I','+'/
      DATA cchar_realtops / ' ',' ','1','1','I','I'/

!     ------ End duplicate definitions common to wrapper routine

       ncolprint=0

      CALL SCOPS(
     &     npoints,
     &     nlev,
     &     ncol,
     &     seed,
     &     cc,
     &     conv,
     &     overlap,
     &     frac_out,
     &     ncolprint
     &)

      CALL ICARUS(
     &     debug,
     &     debugcol,
     &     npoints,
     &     sunlit,
     &     nlev,
     &     ncol,
     &     pfull,
     &     phalf,
     &     qv,
     &     cc,
     &     conv,
     &     dtau_s,
     &     dtau_c,
     &     top_height,
     &     top_height_direction,
     &     overlap,
     &     frac_out,
     &     skt,
     &     emsfc_lw,
     &     at,
     &     dem_s,
     &     dem_c,
     &     fq_isccp,
     &     totalcldarea,
     &     meanptop,
     &     meantaucld,
     &     meanalbedocld,
     &     meantb,
     &     meantbclr,
     &     boxtau,
     &     boxptop
     &)

      return
      end

