! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
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
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Jul 2008 - A. Bodas-Salcedo - Added definitions of ISCCP axes
! Oct 2008 - H. Chepfer       - Added PARASOL_NREFL
! Jun 2010 - R. Marchand      - Modified to support quickbeam V3, added ifdef for  
!                               hydrometeor definitions
! May 2015 - D. Swales        - Tidied up. Set up appropriate fields during initialization. 
! June 2015- D. Swales        - Moved hydrometeor class variables to hydro_class_init in
!                               the module quickbeam_optics.
! Mar 2016 - D. Swales        - Added scops_ccfrac. Was previously hardcoded in prec_scops.f90.  
! Mar 2018 - R. Guzman        - Added LIDAR_NTYPE for the OPAQ diagnostics
! Apr 2018 - R. Guzman        - Added parameters for GROUND LIDAR and ATLID simulators
! Nov 2018 - T. Michibata     - Added CloudSat+MODIS Warmrain Diagnostics
! Mar 2024 - C. Wall          - Added MODIS joint-histogram diagnostics
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE MOD_COSP_CONFIG
    USE COSP_KINDS, ONLY: wp,dp
    IMPLICIT NONE

   ! #####################################################################################
   ! Common COSP information
   ! #####################################################################################
    character(len=32) ::   &
         COSP_VERSION              ! COSP Version ID (set in cosp_interface_init)
    real(wp),parameter ::  &
         R_UNDEF      = -1.0E30, & ! Missing value
         R_GROUND     = -1.0E20, & ! Flag for below ground results
         scops_ccfrac = 0.05       ! Fraction of column (or subcolumn) covered with convective
                                   ! precipitation (default is 5%). *NOTE* This quantity may vary
                                   ! between modeling centers.
    logical :: &
         use_vgrid                 ! True=Use new grid for L3 CLOUDAT and CALIPSO
    integer,parameter ::   &
         SR_BINS = 15,           & ! Number of bins (backscattering coefficient) in CALOPSO LIDAR simulator.
         N_HYDRO = 9               ! Number of hydrometeor classes used by quickbeam radar simulator.

    ! ####################################################################################  
    ! Joint histogram bin-boundaries
    ! tau is used by ISCCP and MISR
    ! pres is used by ISCCP
    ! hgt is used by MISR
    ! ReffLiq is used by MODIS
    ! ReffIce is used by MODIS
    ! *NOTE* ALL JOINT-HISTOGRAM BIN BOUNDARIES ARE DECLARED AND DEFINED HERE IN
    !        COSP_CONFIG, WITH THE EXCEPTION OF THE TAU AXIS USED BY THE MODIS SIMULATOR,
    !        WHICH IS SET DURING INITIALIZATION IN COSP_INTERFACE_INIT.
    ! ####################################################################################
    ! Optical depth bin axis
    integer,parameter :: &
         ntau=7  
    real(wp),parameter,dimension(ntau+1) :: &
       tau_binBounds = (/0.0, 0.3, 1.3, 3.6, 9.4, 23., 60., 10000./)
    real(wp),parameter,dimension(ntau) :: &
         tau_binCenters = (/0.15, 0.80, 2.45, 6.5, 16.2, 41.5, 100.0/)
    real(wp),parameter,dimension(2,ntau) :: &
         tau_binEdges = reshape(source=(/0.0, 0.3,  0.3,  1.3,  1.3,  3.6,      3.6,     &
                                         9.4, 9.4, 23.0, 23.0, 60.0, 60.0, 100000.0/),   &
                                         shape=(/2,ntau/)) 

    ! Optical depth bin axes (ONLY USED BY MODIS SIMULATOR IN v1.4)
    integer :: l,k
    integer,parameter :: &
         ntauV1p4 = 6
    real(wp),parameter,dimension(ntauV1p4+1) :: &
         tau_binBoundsV1p4 = (/0.3, 1.3, 3.6, 9.4, 23., 60., 10000./)
    real(wp),parameter,dimension(2,ntauV1p4) :: &
         tau_binEdgesV1p4 = reshape(source =(/tau_binBoundsV1p4(1),((tau_binBoundsV1p4(k),l=1,2),   &
                                             k=2,ntauV1p4),100000._wp/),shape = (/2,ntauV1p4/)) 
    real(wp),parameter,dimension(ntauV1p4) :: &
         tau_binCentersV1p4 = (tau_binEdgesV1p4(1,:)+tau_binEdgesV1p4(2,:))/2._wp  
    
    ! Cloud-top height pressure bin axis
    integer,parameter :: &
         npres = 7     
    real(wp),parameter,dimension(npres+1) :: &
         pres_binBounds = (/0., 180., 310., 440., 560., 680., 800., 10000./)
    real(wp),parameter,dimension(npres) :: &
         pres_binCenters = (/90000., 74000., 62000., 50000., 37500., 24500., 9000./)   
    real(wp),parameter,dimension(2,npres) :: &
         pres_binEdges = reshape(source=(/100000.0, 80000.0, 80000.0, 68000.0, 68000.0,    &
                                           56000.0, 56000.0, 44000.0, 44000.0, 31000.0,    &
                                           31000.0, 18000.0, 18000.0,     0.0/),           &
                                           shape=(/2,npres/))

    ! Cloud-top height bin axis #1
    integer,parameter :: &
         nhgt = 16
    real(wp),parameter,dimension(nhgt+1) :: &
         hgt_binBounds = (/-.99,0.,0.5,1.,1.5,2.,2.5,3.,4.,5.,7.,9.,11.,13.,15.,17.,99./)
    real(wp),parameter,dimension(nhgt) :: &
         hgt_binCenters = 1000*(/0.,0.25,0.75,1.25,1.75,2.25,2.75,3.5,4.5,6.,8.,10.,12.,   &
         14.5,16.,18./)  
    real(wp),parameter,dimension(2,nhgt) :: &
         hgt_binEdges = 1000.0*reshape(source=(/-99.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.5,  &
                                                  1.5, 2.0, 2.0, 2.5, 2.5, 3.0, 3.0, 4.0,  &
                                                  4.0, 5.0, 5.0, 7.0, 7.0, 9.0, 9.0,11.0,  &
                                                  11.0,13.0,13.0,15.0,15.0,17.0,17.0,99.0/),&
                                                  shape=(/2,nhgt/))    

    ! Liquid and Ice particle bins for MODIS joint histogram of optical-depth and particle
    ! size
    ! Bin edges match Pincus et al. 2023 observational data (doi:10.5194/essd-15-2483-2023)
    integer :: i,j
    integer,parameter :: &
         nReffLiq = 6, & ! Number of ReffLiq bins for tau/ReffLiq and LWP/ReffLiq joint-histogram
         nReffIce = 6    ! Number of ReffIce bins for tau/ReffICE and IWP/ReffIce joint-histogram
    real(wp),parameter,dimension(nReffLiq+1) :: &
         reffLIQ_binBounds = (/4.0e-6, 8e-6, 1.0e-5, 1.25e-5, 1.5e-5, 2.0e-5, 3.0e-5/)
    real(wp),parameter,dimension(nReffIce+1) :: &
         reffICE_binBounds = (/5.0e-6, 1.0e-5, 2.0e-5, 3.0e-5, 4.0e-5, 5.0e-5, 6.0e-5/)
    real(wp),parameter,dimension(2,nReffICE) :: &
         reffICE_binEdges = reshape(source=(/reffICE_binBounds(1),((reffICE_binBounds(k),  &
                                    l=1,2),k=2,nReffICE),reffICE_binBounds(nReffICE+1)/),  &
                                    shape = (/2,nReffICE/)) 
    real(wp),parameter,dimension(2,nReffLIQ) :: &
         reffLIQ_binEdges = reshape(source=(/reffLIQ_binBounds(1),((reffLIQ_binBounds(k),  &
                                    l=1,2),k=2,nReffLIQ),reffLIQ_binBounds(nReffLIQ+1)/),  &
                                    shape = (/2,nReffLIQ/))             
    real(wp),parameter,dimension(nReffICE) :: &
         reffICE_binCenters = (reffICE_binEdges(1,:)+reffICE_binEdges(2,:))/2._wp
    real(wp),parameter,dimension(nReffLIQ) :: &
         reffLIQ_binCenters = (reffLIQ_binEdges(1,:)+reffLIQ_binEdges(2,:))/2._wp

    ! LWP and IWP bins for MODIS joint histogram of (1) LWP and droplet size
    ! and (2) IWP and particle size
    integer, parameter :: &
         nLWP = 7, & ! Number of bins for LWP/ReffLiq joint-histogram
         nIWP = 7    ! Number of bins for IWP/ReffIce joint-histogram
    real(wp),parameter,dimension(nLWP+1) :: &
         LWP_binBounds = (/0., 0.01, 0.03, 0.06, 0.10, 0.15, 0.25, 20.0/) ! kg/m2
    real(wp),parameter,dimension(nIWP+1) :: &
         IWP_binBounds = (/0., 0.02, 0.05, 0.10, 0.20, 0.40, 1.00, 20.0/) ! kg/m2
    real(wp),parameter,dimension(2,nLWP) :: &
         LWP_binEdges = reshape(source=(/LWP_binBounds(1),((LWP_binBounds(k),  &
                                l=1,2),k=2,nLWP),LWP_binBounds(nLWP+1)/),  &
                                shape = (/2,nLWP/))
    real(wp),parameter,dimension(2,nIWP) :: &
         IWP_binEdges = reshape(source=(/IWP_binBounds(1),((IWP_binBounds(k),  &
                                l=1,2),k=2,nIWP),IWP_binBounds(nIWP+1)/),  &
                                shape = (/2,nIWP/))
    real(wp),parameter,dimension(nLWP) :: &
         LWP_binCenters = (LWP_binEdges(1,:)+LWP_binEdges(2,:))/2._wp
    real(wp),parameter,dimension(nIWP) :: &
         IWP_binCenters = (IWP_binEdges(1,:)+IWP_binEdges(2,:))/2._wp
    ! ####################################################################################  
    ! Constants used by RTTOV.
    ! ####################################################################################  
    integer,parameter :: &
         RTTOV_MAX_CHANNELS = 20
    character(len=256),parameter :: &
         rttovDir = '/homedata/rguzman/CALIPSO/RTTOV/rttov_11.3/'
    ! ####################################################################################  
    ! Constants used by the PARASOL simulator.
    ! ####################################################################################  
    integer,parameter :: &
         PARASOL_NREFL = 5,  & ! Number of angles in LUT
         PARASOL_NTAU  = 7     ! Number of optical depths in LUT
    real(wp),parameter,dimension(PARASOL_NREFL) :: &
         PARASOL_SZA = (/0.0, 20.0, 40.0, 60.0, 80.0/)
    REAL(WP),parameter,dimension(PARASOL_NTAU) :: &
         PARASOL_TAU = (/0., 1., 5., 10., 20., 50., 100./)
    
    ! LUTs
    REAL(WP),parameter,dimension(PARASOL_NREFL,PARASOL_NTAU) :: &
         ! LUT for liquid particles
         rlumA = reshape(source=(/ 0.03,     0.03,     0.03,     0.03,     0.03,         &
                                   0.090886, 0.072185, 0.058410, 0.052498, 0.034730,     &
                                   0.283965, 0.252596, 0.224707, 0.175844, 0.064488,     &
                                   0.480587, 0.436401, 0.367451, 0.252916, 0.081667,     &
                                   0.695235, 0.631352, 0.509180, 0.326551, 0.098215,     &
                                   0.908229, 0.823924, 0.648152, 0.398581, 0.114411,     &
                                   1.0,      0.909013, 0.709554, 0.430405, 0.121567/),   &
                                   shape=(/PARASOL_NREFL,PARASOL_NTAU/)),                & 
         ! LUT for ice particles         			     
         rlumB = reshape(source=(/ 0.03,     0.03,     0.03,     0.03,     0.03,         &
                                   0.092170, 0.087082, 0.083325, 0.084935, 0.054157,     &
                                   0.311941, 0.304293, 0.285193, 0.233450, 0.089911,     &
                                   0.511298, 0.490879, 0.430266, 0.312280, 0.107854,     &
                                   0.712079, 0.673565, 0.563747, 0.382376, 0.124127,     &
                                   0.898243, 0.842026, 0.685773, 0.446371, 0.139004,     &
                                   0.976646, 0.912966, 0.737154, 0.473317, 0.145269/),   &
                                   shape=(/PARASOL_NREFL,PARASOL_NTAU/))  

    ! ####################################################################################
    ! ISCCP simulator tau/CTP joint histogram information
    ! ####################################################################################
    integer,parameter :: &
         numISCCPTauBins  = ntau, &              ! Number of optical depth bins
         numISCCPPresBins = npres                ! Number of pressure bins     
    real(wp),parameter,dimension(ntau+1) :: &
         isccp_histTau = tau_binBounds           ! Joint-histogram boundaries (optical depth)
    real(wp),parameter,dimension(npres+1) :: &
         isccp_histPres = pres_binBounds         ! Joint-histogram boundaries (cloud pressure)
    real(wp),parameter,dimension(ntau) :: &
         isccp_histTauCenters = tau_binCenters   ! Joint histogram bin centers (optical depth)
    real(wp),parameter,dimension(npres) :: &   
         isccp_histPresCenters = pres_binCenters ! Joint histogram bin centers (cloud pressure) 
    real(wp),parameter,dimension(2,ntau) :: &
         isccp_histTauEdges = tau_binEdges       ! Joint histogram bin edges (optical depth)
    real(wp),parameter,dimension(2,npres) :: &    
         isccp_histPresEdges = pres_binEdges     ! Joint histogram bin edges (cloud pressure)   
    
    ! ####################################################################################
    ! MISR simulator tau/CTH joint histogram information 
    ! ####################################################################################
    integer,parameter ::  &
         numMISRHgtBins = nhgt, &             ! Number of cloud-top height bins
         numMISRTauBins = ntau                ! Number of optical depth bins
    ! Joint histogram boundaries
    real(wp),parameter,dimension(numMISRHgtBins+1) :: &
         misr_histHgt = hgt_binBounds         ! Joint-histogram boundaries (cloud height)
    real(wp),parameter,dimension(numMISRTauBins+1) :: &
         misr_histTau = tau_binBounds         ! Joint-histogram boundaries (optical-depth)
    real(wp),parameter,dimension(numMISRHgtBins) :: &
         misr_histHgtCenters = hgt_binCenters ! Joint-histogram bin centers (cloud height)
    real(wp),parameter,dimension(2,numMISRHgtBins) :: &
         misr_histHgtEdges = hgt_BinEdges     ! Joint-histogram bin edges (cloud height)
 
    ! ####################################################################################
    ! MODIS simulator tau/CTP joint histogram information 
    ! ####################################################################################
    integer,parameter :: &
         numMODISPresBins = npres                    ! Number of pressure bins for joint-histogram    
    real(wp),parameter,dimension(numMODISPresBins + 1) :: & 
         modis_histPres = 100*pres_binBounds         ! Joint-histogram boundaries (cloud pressure)
    real(wp),parameter,dimension(2, numMODISPresBins) :: &
         modis_histPresEdges = 100*pres_binEdges     ! Joint-histogram bin edges (cloud pressure)
    real(wp),parameter,dimension(numMODISPresBins) :: &
         modis_histPresCenters = 100*pres_binCenters ! Joint-histogram bin centers (cloud pressure)

    ! For the MODIS simulator we want to preserve the ability for cospV1.4.0 to use the
    ! old histogram bin boundaries for optical depth, so these are set up in initialization.
    integer :: &
         numMODISTauBins          ! Number of tau bins for joint-histogram
    real(wp),allocatable,dimension(:) :: &
         modis_histTau            ! Joint-histogram boundaries (optical depth)
    real(wp),allocatable,dimension(:,:) :: &
         modis_histTauEdges       ! Joint-histogram bin edges (optical depth)
    real(wp),allocatable,dimension(:) :: &
         modis_histTauCenters     ! Joint-histogram bin centers (optical depth)
    
    ! ####################################################################################
    ! MODIS simulator tau/ReffICE and tau/ReffLIQ joint-histogram information
    ! ####################################################################################
    ! Ice
    integer,parameter :: &
         numMODISReffIceBins = nReffIce                ! Number of bins for joint-histogram
    real(wp),parameter,dimension(nReffIce+1) :: &
         modis_histReffIce = reffICE_binBounds         ! Effective radius bin boundaries
    real(wp),parameter,dimension(nReffIce) :: &
         modis_histReffIceCenters = reffICE_binCenters ! Effective radius bin centers
    real(wp),parameter,dimension(2,nReffICE) :: &
         modis_histReffIceEdges = reffICE_binEdges     ! Effective radius bin edges
       
    ! Liquid
    integer,parameter :: &
         numMODISReffLiqBins = nReffLiq                ! Number of bins for joint-histogram
    real(wp),parameter,dimension(nReffLiq+1) :: &
         modis_histReffLiq = reffLIQ_binBounds         ! Effective radius bin boundaries 
    real(wp),parameter,dimension(nReffLiq) :: &
         modis_histReffLiqCenters = reffLIQ_binCenters ! Effective radius bin centers
    real(wp),parameter,dimension(2,nReffLiq) :: &
         modis_histReffLiqEdges = reffLIQ_binEdges     ! Effective radius bin edges
    ! ####################################################################################
    ! MODIS simulator LWP/ReffLIQ and IWP/ReffIce joint-histogram information
    ! ####################################################################################
    ! Liquid
    integer,parameter :: &
         numMODISLWPBins = nLWP                        ! Number of bins for joint-histogram
    real(wp),parameter,dimension(nLWP+1) :: &
         modis_histLWP = LWP_binBounds                 ! LWP bin boundaries 
    real(wp),parameter,dimension(nLWP) :: &
         modis_histLWPCenters = LWP_binCenters         ! LWP bin centers
    real(wp),parameter,dimension(2,nLWP) :: &
         modis_histLWPEdges = LWP_binEdges             ! LWP bin edges

    ! Ice
    integer,parameter :: &
         numMODISIWPBins = nIWP                        ! Number of bins for joint-histogram
    real(wp),parameter,dimension(nIWP+1) :: &
         modis_histIWP = IWP_binBounds                 ! IWP bin boundaries 
    real(wp),parameter,dimension(nIWP) :: &
         modis_histIWPCenters = IWP_binCenters         ! IWP bin centers
    real(wp),parameter,dimension(2,nIWP) :: &
         modis_histIWPEdges = IWP_binEdges             ! IWP bin edges     

    ! ####################################################################################
    ! CLOUDSAT reflectivity histogram information 
    ! ####################################################################################
    integer,parameter :: &
       CLOUDSAT_DBZE_BINS     =   15, & ! Number of dBZe bins in histogram (cfad)
       CLOUDSAT_DBZE_MIN      = -100, & ! Minimum value for radar reflectivity
       CLOUDSAT_DBZE_MAX      =   80, & ! Maximum value for radar reflectivity
       CLOUDSAT_CFAD_ZE_MIN   =  -50, & ! Lower value of the first CFAD Ze bin
       CLOUDSAT_CFAD_ZE_WIDTH =    5    ! Bin width (dBZe)

    real(wp),parameter,dimension(CLOUDSAT_DBZE_BINS+1) :: &
         cloudsat_histRef = (/CLOUDSAT_DBZE_MIN,(/(i, i=int(CLOUDSAT_CFAD_ZE_MIN+CLOUDSAT_CFAD_ZE_WIDTH),&
                             int(CLOUDSAT_CFAD_ZE_MIN+(CLOUDSAT_DBZE_BINS-1)*CLOUDSAT_CFAD_ZE_WIDTH),    &
                             int(CLOUDSAT_CFAD_ZE_WIDTH))/),CLOUDSAT_DBZE_MAX/)
    real(wp),parameter,dimension(2,CLOUDSAT_DBZE_BINS) :: &
         cloudsat_binEdges = reshape(source=(/cloudsat_histRef(1),((cloudsat_histRef(k), &
                                   l=1,2),k=2,CLOUDSAT_DBZE_BINS),cloudsat_histRef(CLOUDSAT_DBZE_BINS+1)/),&
                                   shape = (/2,CLOUDSAT_DBZE_BINS/))     
    real(wp),parameter,dimension(CLOUDSAT_DBZE_BINS) :: &
         cloudsat_binCenters = (cloudsat_binEdges(1,:)+cloudsat_binEdges(2,:))/2._wp

    ! Parameters for Cloudsat near-surface precipitation diagnostics.
    ! Precipitation classes.
    integer, parameter :: &
         nCloudsatPrecipClass = 10
    integer, parameter :: &
         pClass_noPrecip      = 0, & ! No precipitation
         pClass_Rain1         = 1, & ! Rain possible
         pClass_Rain2         = 2, & ! Rain probable
         pClass_Rain3         = 3, & ! Rain certain
         pClass_Snow1         = 4, & ! Snow possible
         pClass_Snow2         = 5, & ! Snow certain
         pClass_Mixed1        = 6, & ! Mixed-precipitation possible
         pClass_Mixed2        = 7, & ! Mixed-precipitation certain
         pClass_Rain4         = 8, & ! Heavy rain
         pClass_default       = 9    ! Default
    ! Reflectivity bin boundaries, used by decision tree to classify precipitation type.
    real(wp), dimension(4),parameter :: &
         Zenonbinval =(/0._wp, -5._wp, -7.5_wp, -15._wp/)
    real(wp), dimension(6),parameter :: &
         Zbinvallnd = (/10._wp, 5._wp, 2.5_wp, -2.5_wp, -5._wp, -15._wp/)
    ! Vertical level index(Nlvgrid) for Cloudsat precipitation occurence/frequency diagnostics.
    ! Level 39 of Nlvgrid(40) is 480-960m.
    integer, parameter :: &
         cloudsat_preclvl = 39

    ! ####################################################################################
    ! CLOUDSAT and MODIS joint product information (2018.11.22)
    ! ####################################################################################
    ! @ COSP_DIAG_WARMRAIN:
    integer, parameter :: CFODD_NCLASS  =    3 ! # of classes for CFODD (classified by MODIS Reff)
    integer, parameter :: WR_NREGIME    =    3 ! # of warm-rain regimes (non-precip/drizzling/raining)
    integer, parameter :: SGCLD_CLR     =    0 ! sub-grid cloud ID (fracout): clear-sky
    integer, parameter :: SGCLD_ST      =    1 ! sub-grid cloud ID (fracout): stratiform
    integer, parameter :: SGCLD_CUM     =    2 ! sub-grid cloud ID (fracout): cumulus
    real(wp),parameter :: CWP_THRESHOLD = 0.00 ! cloud water path threshold
    real(wp),parameter :: COT_THRESHOLD = 0.30 ! cloud optical thickness threshold
    real(wp),parameter,dimension(CFODD_NCLASS+1) :: &
         CFODD_BNDRE = (/5.0e-6, 12.0e-6, 18.0e-6, 35.0e-6/) ! Reff bnds
    real(wp),parameter,dimension(2) :: &
         CFODD_BNDZE = (/-15.0, 0.0/)                        ! dBZe bnds (cloud/drizzle/precip)
    real(wp),parameter :: CFODD_DBZE_MIN    =  -30.0 ! Minimum value of CFODD dBZe bin
    real(wp),parameter :: CFODD_DBZE_MAX    =   20.0 ! Maximum value of CFODD dBZe bin
    real(wp),parameter :: CFODD_ICOD_MIN    =    0.0 ! Minimum value of CFODD ICOD bin
    real(wp),parameter :: CFODD_ICOD_MAX    =   60.0 ! Maximum value of CFODD ICOD bin
    real(wp),parameter :: CFODD_DBZE_WIDTH  =    2.0 ! Bin width (dBZe)
    real(wp),parameter :: CFODD_ICOD_WIDTH  =    2.0 ! Bin width (ICOD)
    integer,parameter :: &
         CFODD_NDBZE = INT( (CFODD_DBZE_MAX-CFODD_DBZE_MIN)/CFODD_DBZE_WIDTH ) ! Number of CFODD dBZe bins
    integer,parameter :: &
         CFODD_NICOD = INT( (CFODD_ICOD_MAX-CFODD_ICOD_MIN)/CFODD_ICOD_WIDTH ) ! Number of CFODD ICOD bins
    real(wp),parameter,dimension(CFODD_NDBZE+1) :: &
         CFODD_HISTDBZE = (/int(CFODD_DBZE_MIN),(/(i, i=int(CFODD_DBZE_MIN+CFODD_DBZE_WIDTH), &
                           int(CFODD_DBZE_MIN+(CFODD_NDBZE-1)*CFODD_DBZE_WIDTH),              &
                           int(CFODD_DBZE_WIDTH))/),int(CFODD_DBZE_MAX)/)
    real(wp),parameter,dimension(CFODD_NICOD+1) :: &
         CFODD_HISTICOD = (/int(CFODD_ICOD_MIN),(/(i, i=int(CFODD_ICOD_MIN+CFODD_ICOD_WIDTH), &
                           int(CFODD_ICOD_MIN+(CFODD_NICOD-1)*CFODD_ICOD_WIDTH),              &
                           int(CFODD_ICOD_WIDTH))/),int(CFODD_ICOD_MAX)/)
    real(wp),parameter,dimension(2,CFODD_NDBZE) :: &
         CFODD_HISTDBZEedges = reshape(source=(/CFODD_HISTDBZE(1),((CFODD_HISTDBZE(k),    &
                                 l=1,2),k=2,CFODD_NDBZE),CFODD_HISTDBZE(CFODD_NDBZE+1)/), &
                                 shape = (/2,CFODD_NDBZE/))
    real(wp),parameter,dimension(CFODD_NDBZE) :: &
         CFODD_HISTDBZEcenters = (CFODD_HISTDBZEedges(1,:)+CFODD_HISTDBZEedges(2,:))/2._wp
    real(wp),parameter,dimension(2,CFODD_NICOD) :: &
         CFODD_HISTICODedges = reshape(source=(/CFODD_HISTICOD(1),((CFODD_HISTICOD(k),    &
                                 l=1,2),k=2,CFODD_NICOD),CFODD_HISTICOD(CFODD_NICOD+1)/), &
                                 shape = (/2,CFODD_NICOD/))
    real(wp),parameter,dimension(CFODD_NICOD) :: &
         CFODD_HISTICODcenters = (CFODD_HISTICODedges(1,:)+CFODD_HISTICODedges(2,:))/2._wp

    ! ####################################################################################
    ! Parameters used by the CALIPSO LIDAR simulator
    ! #################################################################################### 
    ! CALISPO backscatter histogram bins 
    real(wp),parameter ::     &
       S_cld       = 5.0,     & ! Threshold for cloud detection
       S_att       = 0.01,    & !
       S_cld_att   = 30.        ! Threshold for undefined cloud phase detection
    real(wp),parameter,dimension(SR_BINS+1) :: &
         calipso_histBsct = (/-1.,0.01,1.2,3.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,   &
                              60.0,80.0,999./)         ! Backscatter histogram bins
    real(wp),parameter,dimension(2,SR_BINS) :: &
         calipso_binEdges = reshape(source=(/calipso_histBsct(1),((calipso_histBsct(k),  &
                                    l=1,2),k=2,SR_BINS),calipso_histBsct(SR_BINS+1)/),   &
                                    shape = (/2,SR_BINS/))     
    real(wp),parameter,dimension(SR_BINS) :: &
         calipso_binCenters = (calipso_binEdges(1,:)+calipso_binEdges(2,:))/2._wp  

    integer,parameter  ::     &
       LIDAR_NTEMP = 40, & 
       LIDAR_NCAT  = 4,  & ! Number of categories for cloudtop heights (high/mid/low/tot)
       LIDAR_NTYPE = 3     ! Number of categories for OPAQ (opaque/thin cloud + z_opaque)
    real(wp),parameter,dimension(LIDAR_NTEMP) :: &
       LIDAR_PHASE_TEMP=                                                                 &
       (/-91.5,-88.5,-85.5,-82.5,-79.5,-76.5,-73.5,-70.5,-67.5,-64.5,                    &
         -61.5,-58.5,-55.5,-52.5,-49.5,-46.5,-43.5,-40.5,-37.5,-34.5,                    &
         -31.5,-28.5,-25.5,-22.5,-19.5,-16.5,-13.5,-10.5, -7.5, -4.5,                    &
          -1.5,  1.5,  4.5,  7.5, 10.5, 13.5, 16.5, 19.5, 22.5, 25.5/)
    real(wp),parameter,dimension(2,LIDAR_NTEMP) :: &
       LIDAR_PHASE_TEMP_BNDS=reshape(source=                                             &
          (/-273.15, -90., -90., -87., -87., -84., -84., -81., -81., -78.,               &
             -78.,   -75., -75., -72., -72., -69., -69., -66., -66., -63.,               &
             -63.,   -60., -60., -57., -57., -54., -54., -51., -51., -48.,               &
             -48.,   -45., -45., -42., -42., -39., -39., -36., -36., -33.,               &
             -33.,   -30., -30., -27., -27., -24., -24., -21., -21., -18.,               &
             -18.,   -15., -15., -12., -12.,  -9.,  -9.,  -6.,  -6.,  -3.,               &
              -3.,     0.,   0.,   3.,   3.,   6.,   6.,   9.,   9.,  12.,               &
              12.,    15.,  15.,  18.,  18.,  21.,  21.,  24.,  24., 100. /),            &
              shape=(/2,40/))        

    ! ####################################################################################
    ! Parameters used by the GROUND LIDAR simulator
    ! #################################################################################### 
    ! GROUND LIDAR backscatter histogram bins 
!    real(wp),parameter ::     &
!       S_cld       = 5.0,     & ! Threshold for cloud detection
!       S_att       = 0.01,    & !
!       S_cld_att   = 30.        ! Threshold for undefined cloud phase detection
    real(wp),parameter,dimension(SR_BINS+1) :: &
         grLidar532_histBsct = (/-1.,0.01,1.2,3.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,  &
                                 60.0,80.0,999./)         ! Backscatter histogram bins
    real(wp),parameter,dimension(2,SR_BINS) :: &
         grLidar532_binEdges = reshape(source=(/grLidar532_histBsct(1),((grLidar532_histBsct(k),  &
                                    l=1,2),k=2,SR_BINS),grLidar532_histBsct(SR_BINS+1)/),   &
                                    shape = (/2,SR_BINS/))     
    real(wp),parameter,dimension(SR_BINS) :: &
         grLidar532_binCenters = (grLidar532_binEdges(1,:)+grLidar532_binEdges(2,:))/2._wp  

!    integer,parameter  ::     &
!       LIDAR_NCAT  = 4       ! Number of categories for cloudtop heights (high/mid/low/tot)

    ! ####################################################################################
    ! Parameters used by the ATLID LIDAR simulator
    ! #################################################################################### 
    ! ATLID LIDAR backscatter histogram bins 
    real(wp),parameter ::     &
       S_cld_atlid       = 1.74,    & ! Threshold for cloud detection
       S_att_atlid       = 0.01,    & !
       S_cld_att_atlid   = 6.67        ! Threshold for undefined cloud phase detection
    real(wp),parameter,dimension(SR_BINS+1) :: &
         atlid_histBsct = (/-1.,0.01,1.03,1.38,1.74,2.07,2.62,3.65,4.63,5.63,6.67,8.8,11.25,  &
                                 13.2,17.2,999./)         ! Backscatter histogram bins
    real(wp),parameter,dimension(2,SR_BINS) :: &
         atlid_binEdges = reshape(source=(/atlid_histBsct(1),((atlid_histBsct(k),  &
                                    l=1,2),k=2,SR_BINS),atlid_histBsct(SR_BINS+1)/),   &
                                    shape = (/2,SR_BINS/))     
    real(wp),parameter,dimension(SR_BINS) :: &
         atlid_binCenters = (atlid_binEdges(1,:)+atlid_binEdges(2,:))/2._wp  

!    integer,parameter  ::     &
!       LIDAR_NCAT  = 4       ! Number of categories for cloudtop heights (high/mid/low/tot)

    ! ####################################################################################
    ! New vertical grid used by CALIPSO and CLOUDSAT L3 (set up during initialization)
    ! ####################################################################################
    integer :: &
         Nlvgrid      ! Number of levels in New grid
    real(wp),dimension(:),allocatable :: &
       vgrid_zl,  & ! New grid bottoms
       vgrid_zu,  & ! New grid tops
       vgrid_z,   & ! New grid center
       dz           ! dZ

END MODULE MOD_COSP_CONFIG
