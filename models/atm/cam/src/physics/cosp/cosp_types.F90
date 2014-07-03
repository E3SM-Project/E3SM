! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Feb 2008 - R. Marchand      - Added Quickbeam types and initialisation
! Oct 2008 - H. Chepfer       - Added PARASOL reflectance diagnostic
! Nov 2008 - R. Marchand      - Added MISR diagnostics
! Nov 2008 - V. John          - Added RTTOV diagnostics
!
! 
MODULE MOD_COSP_TYPES
    USE MOD_COSP_CONSTANTS
    USE MOD_COSP_UTILS

    use radar_simulator_types, only: class_param, mie, nd, mt_nd, dmax, dmin, mt_ttl, mt_tti, cnt_liq, cnt_ice  ! added by roj Feb 2008

    IMPLICIT NONE
    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------- DERIVED TYPES ----------------------------    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Configuration choices (simulators, variables)
  TYPE COSP_CONFIG
     logical :: Lradar_sim,Llidar_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim,Lstats,Lwrite_output, &
                Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
                LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp,Lcllcalipso, &
                Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lcltradar,Lcltradar2,Lpctisccp,Ldbze94,Ltauisccp,Lcltisccp, & !modified by ZYY
                Llongitude,Llatitude,LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
                Lfracout,LlidarBetaMol532,Ltbrttov, &
                Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
                Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
                Liwpmodis,Lclmodis

     character(len=32) :: out_list(N_OUT_LIST)
  END TYPE COSP_CONFIG
  
  ! Outputs from RTTOV
  TYPE COSP_RTTOV
     ! Dimensions
     integer :: Npoints   ! Number of gridpoints
     integer :: Nchan     ! Number of channels
     
     ! Brightness temperatures (Npoints,Nchan)
     real,pointer :: tbs(:,:)
     
  END TYPE COSP_RTTOV
  
  ! Outputs from MISR simulator
  TYPE COSP_MISR
     ! Dimensions
     integer :: Npoints   ! Number of gridpoints
     integer :: Ntau      ! Number of tau intervals
     integer :: Nlevels   ! Number of cth levels

     ! --- (npoints,ntau,nlevels)
     !  the fraction of the model grid box covered by each of the MISR cloud types
     real,pointer :: fq_MISR(:,:,:)  
     
     ! --- (npoints)
     real,pointer :: MISR_meanztop(:), MISR_cldarea(:)
     ! --- (npoints,nlevels)
     real,pointer :: MISR_dist_model_layertops(:,:)
  END TYPE COSP_MISR

  ! Outputs from ISCCP simulator
  TYPE COSP_ISCCP
     ! Dimensions
     integer :: Npoints   ! Number of gridpoints
     integer :: Ncolumns  ! Number of columns
     integer :: Nlevels   ! Number of levels

    
     ! --- (npoints,tau=7,pressure=7)
     !  the fraction of the model grid box covered by each of the 49 ISCCP D level cloud types
     real,pointer :: fq_isccp(:,:,:)
     
     ! --- (npoints) ---
     ! The fraction of model grid box columns with cloud somewhere in them.
     ! This should equal the sum over all entries of fq_isccp
     real,pointer :: totalcldarea(:)
     ! mean all-sky 10.5 micron brightness temperature
     real,pointer ::  meantb(:)
     ! mean clear-sky 10.5 micron brightness temperature
     real,pointer ::  meantbclr(:)
     
     ! The following three means are averages over the cloudy areas only.  If no
     ! clouds are in grid box all three quantities should equal zero.
     
     !  mean cloud top pressure (mb) - linear averaging in cloud top pressure.
     real,pointer :: meanptop(:)
     !  mean optical thickness linear averaging in albedo performed.
     real,pointer :: meantaucld(:)
     ! mean cloud albedo. linear averaging in albedo performed 
     real,pointer :: meanalbedocld(:)  
     
     !--- (npoints,ncol) ---
     !  optical thickness in each column     
     real,pointer :: boxtau(:,:)
     !  cloud top pressure (mb) in each column
     real,pointer :: boxptop(:,:)        
  END TYPE COSP_ISCCP
  
  ! Summary statistics from radar
  TYPE COSP_VGRID
    logical :: use_vgrid ! Logical flag that indicates change of grid
    logical :: csat_vgrid ! Flag for Cloudsat grid
    integer :: Npoints   ! Number of sampled points
    integer :: Ncolumns  ! Number of subgrid columns
    integer :: Nlevels   ! Number of model levels
    integer :: Nlvgrid   ! Number of levels of new grid
    ! Array with dimensions (Nlvgrid)
    real, dimension(:), pointer :: z,zl,zu ! Height and lower and upper boundaries of new levels
    ! Array with dimensions (Nlevels)
    real, dimension(:), pointer :: mz,mzl,mzu ! Height and lower and upper boundaries of model levels
  END TYPE COSP_VGRID
  
  ! Output data from lidar code
  TYPE COSP_SGLIDAR
    ! Dimensions
    integer :: Npoints   ! Number of gridpoints
    integer :: Ncolumns  ! Number of columns
    integer :: Nlevels   ! Number of levels
    integer :: Nhydro    ! Number of hydrometeors    
    integer :: Nrefl     ! Number of parasol reflectances
    ! Arrays with dimensions (Npoints,Nlevels)
    real,dimension(:,:),pointer :: beta_mol   ! Molecular backscatter
    ! Arrays with dimensions (Npoints,Ncolumns,Nlevels)
    real,dimension(:,:,:),pointer :: beta_tot   ! Total backscattered signal
    real,dimension(:,:,:),pointer :: tau_tot    ! Optical thickness integrated from top to level z
    ! Arrays with dimensions (Npoints,Ncolumns,Nrefl)
    real,dimension(:,:,:),pointer :: refl       ! parasol reflectances
  END TYPE COSP_SGLIDAR
  
  ! Output data from radar code
  TYPE COSP_SGRADAR
    ! Dimensions
    integer :: Npoints   ! Number of gridpoints
    integer :: Ncolumns  ! Number of columns
    integer :: Nlevels   ! Number of levels
    integer :: Nhydro    ! Number of hydrometeors
    ! output vertical levels: spaceborne radar -> from TOA to SURFACE
    ! Arrays with dimensions (Npoints,Nlevels)
    real,dimension(:,:),pointer :: att_gas ! 2-way attenuation by gases [dBZ]
    ! Arrays with dimensions (Npoints,Ncolumns,Nlevels)
    real,dimension(:,:,:),pointer :: Ze_tot ! Effective reflectivity factor [dBZ]
 
  END TYPE COSP_SGRADAR

  
  ! Summary statistics from radar
  TYPE COSP_RADARSTATS
    integer :: Npoints  ! Number of sampled points
    integer :: Ncolumns ! Number of subgrid columns
    integer :: Nlevels  ! Number of model levels
    integer :: Nhydro   ! Number of hydrometeors
    ! Array with dimensions (Npoints,dBZe_bins,Nlevels)
    real, dimension(:,:,:), pointer :: cfad_ze ! Ze CFAD
    ! Array with dimensions (Npoints)
    real,dimension(:),pointer :: radar_lidar_tcc ! Radar&lidar total cloud amount, grid-box scale
    real,dimension(:),pointer :: radar_tcc !Radar total cloud amount !modified by ZYY
    real,dimension(:),pointer :: radar_tcc_2 !Radar total cloud amount without the data for the first kilometer above surface !modified by ZYY
    ! Arrays with dimensions (Npoints,Nlevels)
    real, dimension(:,:),pointer :: lidar_only_freq_cloud
  END TYPE COSP_RADARSTATS

  ! Summary statistics from lidar
  TYPE COSP_LIDARSTATS
    integer :: Npoints  ! Number of sampled points
    integer :: Ncolumns ! Number of subgrid columns
    integer :: Nlevels  ! Number of model levels
    integer :: Nhydro   ! Number of hydrometeors
    integer :: Nrefl    ! Number of parasol reflectances
    
    ! Arrays with dimensions (SR_BINS)
    real, dimension(:),pointer :: srbval ! SR bins in cfad_sr
    ! Arrays with dimensions (Npoints,SR_BINS,Nlevels)
    real, dimension(:,:,:),pointer :: cfad_sr   ! CFAD of scattering ratio
    ! Arrays with dimensions (Npoints,Nlevels)
    real, dimension(:,:),pointer :: lidarcld    ! 3D "lidar" cloud fraction 
    ! Arrays with dimensions (Npoints,LIDAR_NCAT)
    real, dimension(:,:),pointer :: cldlayer      ! low, mid, high-level lidar cloud cover
    ! Arrays with dimensions (Npoints,PARASOL_NREFL)
    real, dimension(:,:),pointer :: parasolrefl   ! mean parasol reflectance

  END TYPE COSP_LIDARSTATS

    
  ! Input data for simulator. Subgrid scale.
  ! Input data from SURFACE to TOA
  TYPE COSP_SUBGRID
    ! Dimensions
    integer :: Npoints   ! Number of gridpoints
    integer :: Ncolumns  ! Number of columns
    integer :: Nlevels   ! Number of levels
    integer :: Nhydro    ! Number of hydrometeors
    
    real,dimension(:,:,:),pointer :: prec_frac  ! Subgrid precip array. Dimensions (Npoints,Ncolumns,Nlevels)
    real,dimension(:,:,:),pointer :: frac_out  ! Subgrid cloud array. Dimensions (Npoints,Ncolumns,Nlevels)
  END TYPE COSP_SUBGRID

  ! Input data for simulator at Subgrid scale.
  ! Used on a reduced number of points
  TYPE COSP_SGHYDRO
    ! Dimensions
    integer :: Npoints   ! Number of gridpoints
    integer :: Ncolumns  ! Number of columns
    integer :: Nlevels   ! Number of levels
    integer :: Nhydro    ! Number of hydrometeors
    real,dimension(:,:,:,:),pointer :: mr_hydro ! Mixing ratio of each hydrometeor 
                                                ! (Npoints,Ncolumns,Nlevels,Nhydro) [kg/kg]
    real,dimension(:,:,:,:),pointer :: Reff     ! Effective Radius of each hydrometeor
                                                ! (Reff==0 means use default size)   
                                                ! (Npoints,Ncolumns,Nlevels,Nhydro) [m]
  END TYPE COSP_SGHYDRO
  
  ! Input data for simulator. Gridbox scale.
  TYPE COSP_GRIDBOX
    ! Scalars and dimensions
    integer :: Npoints   ! Number of gridpoints
    integer :: Nlevels   ! Number of levels
    integer :: Ncolumns  ! Number of columns
    integer :: Nhydro    ! Number of hydrometeors
    integer :: Nprmts_max_hydro    ! Max number of parameters for hydrometeor size distributions
    integer :: Naero    ! Number of aerosol species
    integer :: Nprmts_max_aero    ! Max number of parameters for aerosol size distributions
    integer :: Npoints_it   ! Max number of gridpoints to be processed in one iteration
    
    ! Time [days]
    ! jsb double precision -> real*8
    real*8 :: time
    real*8 :: time_bnds(2)
    ! jsb
    
    ! Radar ancillary info
    real :: radar_freq, & ! Radar frequency [GHz]
            k2 ! |K|^2, -1=use frequency dependent default
    integer :: surface_radar, & ! surface=1, spaceborne=0
           use_mie_tables, & ! use a precomputed loopup table? yes=1,no=0
           use_gas_abs, & ! include gaseous absorption? yes=1,no=0
           do_ray, & ! calculate/output Rayleigh refl=1, not=0
           melt_lay ! melting layer model off=0, on=1
 
    ! structures used by radar simulator that need to be set only ONCE per radar configuration (e.g. freq, pointing direction) ... added by roj Feb 2008
    type(class_param) ::  hp    ! structure used by radar simulator to store Ze and N scaling constants and other information
    type(mie)::  mt     ! structure used by radar simulator to store mie LUT information
    integer :: nsizes       ! number of discrete drop sizes (um) used to represent the distribution
    real*8, dimension(:), pointer :: D ! array of discrete drop sizes (um) used to represent the distribution
    real*8, dimension(:), pointer :: mt_ttl, mt_tti ! array of temperatures used with Ze_scaling (also build into mie LUT)
    
    ! Lidar
    integer :: lidar_ice_type !ice particle shape hypothesis in lidar calculations 
                              !(ice_type=0 for spheres, ice_type=1 for non spherical particles)
    
    ! Radar
    logical ::  use_precipitation_fluxes  ! True if precipitation fluxes are input to the algorithm 
    logical ::  use_reff  ! True if Reff is to be used by radar 
    
    ! Geolocation (Npoints)
    real,dimension(:),pointer :: longitude ! longitude [degrees East]
    real,dimension(:),pointer :: latitude  ! latitude [deg North]
    ! Gridbox information (Npoints,Nlevels)
    real,dimension(:,:),pointer :: zlev ! Height of model levels [m]
    real,dimension(:,:),pointer :: zlev_half ! Height at half model levels [m] (Bottom of model layer)
    real,dimension(:,:),pointer :: dlev ! Depth of model levels  [m]
    real,dimension(:,:),pointer :: p  ! Pressure at full model levels [Pa]
    real,dimension(:,:),pointer :: ph ! Pressure at half model levels [Pa]
    real,dimension(:,:),pointer :: T ! Temperature at model levels [K]
    real,dimension(:,:),pointer :: q  ! Relative humidity to water (%)
    real,dimension(:,:),pointer :: sh ! Specific humidity to water [kg/kg]
    real,dimension(:,:),pointer :: dtau_s ! mean 0.67 micron optical depth of stratiform
                                          !  clouds in each model level
                                          !  NOTE:  this the cloud optical depth of only the
                                          !  cloudy part of the grid box, it is not weighted
                                          !  with the 0 cloud optical depth of the clear
                                          !         part of the grid box
    real,dimension(:,:),pointer :: dtau_c !  mean 0.67 micron optical depth of convective
                                          !  clouds in each model level.  Same note applies as in dtau_s.
    ! jsb added from steve's
    real,dimension(:,:),pointer :: dtau_s_snow !  mean 0.67 micron optical depth of stratiform snow (in-snow)
    !jsb
    real,dimension(:,:),pointer :: dem_s  !  10.5 micron longwave emissivity of stratiform
                                          !  clouds in each model level.  Same note applies as in dtau_s.
    real,dimension(:,:),pointer :: dem_c  !  10.5 micron longwave emissivity of convective
                                          !  clouds in each model level.  Same note applies as in dtau_s.
! jsb added from steve's
    real,dimension(:,:),pointer :: dem_s_snow  !  10.5 micron longwave emissivity of stratiform snow (in-snow)
!jsb
    real,dimension(:,:),pointer :: mr_ozone !  Ozone mass mixing ratio [kg/kg]

    ! Point information (Npoints)
    real,dimension(:),pointer :: land !Landmask [0 - Ocean, 1 - Land]
    real,dimension(:),pointer :: psfc !Surface pressure [Pa]
    real,dimension(:),pointer :: sunlit ! (npoints) 1 for day points, 0 for nightime
    real,dimension(:),pointer :: skt  ! Skin temperature (K)
    real,dimension(:),pointer :: sfc_height  ! Surface height [m]
    real,dimension(:),pointer :: u_wind  ! eastward wind [m s-1]
    real,dimension(:),pointer :: v_wind  ! northward wind [m s-1]

    ! TOTAL and CONV cloud fraction for SCOPS
    real,dimension(:,:),pointer :: tca ! Total cloud fraction
    real,dimension(:,:),pointer :: cca ! Convective cloud fraction
    ! Precipitation fluxes on model levels
    real,dimension(:,:),pointer :: rain_ls ! large-scale precipitation flux of rain [kg/m2.s]
    real,dimension(:,:),pointer :: rain_cv ! convective precipitation flux of rain [kg/m2.s]
    real,dimension(:,:),pointer :: snow_ls ! large-scale precipitation flux of snow [kg/m2.s]
    real,dimension(:,:),pointer :: snow_cv ! convective precipitation flux of snow [kg/m2.s]
    real,dimension(:,:),pointer :: grpl_ls ! large-scale precipitation flux of graupel [kg/m2.s]
    ! Hydrometeors concentration and distribution parameters
!     real,dimension(:,:,:),pointer :: fr_hydro ! Fraction of the gridbox occupied by each hydrometeor (Npoints,Nlevels,Nhydro)
    real,dimension(:,:,:),pointer :: mr_hydro ! Mixing ratio of each hydrometeor (Npoints,Nlevels,Nhydro) [kg/kg]
    real,dimension(:,:),pointer   :: dist_prmts_hydro !Distributional parameters for hydrometeors (Nprmts_max_hydro,Nhydro)
    ! Effective radius [m]. (Npoints,Nlevels,Nhydro)
    real,dimension(:,:,:),pointer :: Reff
    ! Aerosols concentration and distribution parameters
    real,dimension(:,:,:),pointer :: conc_aero ! Aerosol concentration for each species (Npoints,Nlevels,Naero)
    integer,dimension(:),pointer :: dist_type_aero ! Particle size distribution type for each aerosol species (Naero)
    real,dimension(:,:,:,:),pointer :: dist_prmts_aero ! Distributional parameters for aerosols 
                                                       ! (Npoints,Nlevels,Nprmts_max_aero,Naero)
    ! ISCCP simulator inputs
    integer :: isccp_top_height !  1 = adjust top height using both a computed
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
    integer :: isccp_top_height_direction ! direction for finding atmosphere pressure level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 ! 1 = find the *lowest* altitude (highest pressure) level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 ! 2 = find the *highest* altitude (lowest pressure) level
                                 ! with interpolated temperature equal to the radiance 
                                 ! determined cloud-top temperature
                                 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                 ! 1 = default setting, and matches all versions of 
                                 ! ISCCP simulator with versions numbers 3.5.1 and lower
                                 ! 2 = experimental setting  
    integer :: isccp_overlap !  overlap type (1=max, 2=rand, 3=max/rand)
    real :: isccp_emsfc_lw      ! 10.5 micron emissivity of surface (fraction)
  
    ! RTTOV inputs/options
    integer :: plat      ! satellite platform
    integer :: sat       ! satellite
    integer :: inst      ! instrument
    integer :: Nchan     ! Number of channels to be computed
    integer, dimension(:), pointer :: Ichan   ! Channel numbers
    real,    dimension(:), pointer :: Surfem  ! Surface emissivity
    real    :: ZenAng ! Satellite Zenith Angles
    real :: co2,ch4,n2o,co ! Mixing ratios of trace gases

  END TYPE COSP_GRIDBOX
 
CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_RTTOV -------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_RTTOV(Npoints,Nchan,x)
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Nchan ! Number of channels
    type(cosp_rttov),intent(out) :: x
    
    ! Dimensions
    x%Npoints  = Npoints
    x%Nchan    = Nchan
      
    ! --- Allocate arrays ---
    allocate(x%tbs(Npoints, Nchan))
    ! --- Initialise to zero ---
    x%tbs     = 0.0
  END SUBROUTINE CONSTRUCT_COSP_RTTOV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_RTTOV ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_RTTOV(x)
    type(cosp_rttov),intent(inout) :: x
    
    ! --- Deallocate arrays ---
    deallocate(x%tbs)
  END SUBROUTINE FREE_COSP_RTTOV
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_MISR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_MISR(cfg,Npoints,x)
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints   ! Number of gridpoints
    type(cosp_misr),intent(out) :: x
    ! Local variables
    integer :: i,j,k
    
   
    ! Allocate minumum storage if simulator not used
    if (cfg%Lmisr_sim) then
      i = Npoints
      j = 7
      k = MISR_N_CTH
    else
      i = 1
      j = 1
      k = 1
    endif
    
    ! Dimensions
    x%Npoints = i
    x%Ntau    = j
    x%Nlevels = k
    
    ! allocate space for MISR simulator outputs ...
    allocate(x%fq_MISR(i,j,k), x%MISR_meanztop(i),x%MISR_cldarea(i), x%MISR_dist_model_layertops(i,k))
    x%fq_MISR = 0.0
    x%MISR_meanztop = 0.0
    x%MISR_cldarea = 0.0
    x%MISR_dist_model_layertops = 0.0
    
  END SUBROUTINE CONSTRUCT_COSP_MISR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_MISR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_MISR(x)
    type(cosp_misr),intent(inout) :: x
    deallocate(x%fq_MISR, x%MISR_meanztop,x%MISR_cldarea, x%MISR_dist_model_layertops)
    
  END SUBROUTINE FREE_COSP_MISR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_ISCCP ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_ISCCP(cfg,Npoints,Ncolumns,Nlevels,x)
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Ncolumns ! Number of subgrid columns
    integer,intent(in) :: Nlevels  ! Number of model levels
    type(cosp_isccp),intent(out) :: x
    ! Local variables
    integer :: i,j,k
    
    ! Allocate minumum storage if simulator not used
    if (cfg%Lisccp_sim) then
      i = Npoints
      j = Ncolumns
      k = Nlevels
    else
      i = 1
      j = 1
      k = 1
    endif
    
    ! Dimensions
    x%Npoints  = i
    x%Ncolumns = j
    x%Nlevels  = k
    
    ! --- Allocate arrays ---
    allocate(x%fq_isccp(i,7,7), x%totalcldarea(i), &
         x%meanptop(i), x%meantaucld(i), &
         x%meantb(i), x%meantbclr(i), &
         x%boxtau(i,j), x%boxptop(i,j), &
         x%meanalbedocld(i))
    ! --- Initialise to zero ---
    x%fq_isccp     = 0.0
    x%totalcldarea = 0.0
    x%meanptop     = 0.0
    x%meantaucld   = 0.0
    x%meantb       = 0.0
    x%meantbclr    = 0.0
    x%boxtau       = 0.0
    x%boxptop      = 0.0
    x%meanalbedocld= 0.0
  END SUBROUTINE CONSTRUCT_COSP_ISCCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_ISCCP -----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_ISCCP(x)
    type(cosp_isccp),intent(inout) :: x
    
    deallocate(x%fq_isccp, x%totalcldarea, &
         x%meanptop, x%meantaucld, x%meantb, x%meantbclr, &
         x%boxtau, x%boxptop, x%meanalbedocld)
  END SUBROUTINE FREE_COSP_ISCCP
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_VGRID ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_VGRID(gbx,Nlvgrid,use_vgrid,cloudsat,x)
    type(cosp_gridbox),intent(in) :: gbx ! Gridbox information
    integer,intent(in) :: Nlvgrid  ! Number of new levels    
    logical,intent(in) :: use_vgrid! Logical flag that controls the output on a different grid
    logical,intent(in) :: cloudsat ! TRUE if a CloudSat like grid (480m) is requested
    type(cosp_vgrid),intent(out) :: x
    
    ! Local variables
    integer :: i
    real :: zstep
    
    x%use_vgrid  = use_vgrid
    x%csat_vgrid = cloudsat
    
    ! Dimensions
    x%Npoints  = gbx%Npoints
    x%Ncolumns = gbx%Ncolumns
    x%Nlevels  = gbx%Nlevels
    
    ! --- Allocate arrays ---
    if (use_vgrid) then
      x%Nlvgrid = Nlvgrid
    else 
      x%Nlvgrid = gbx%Nlevels
    endif
    allocate(x%z(x%Nlvgrid),x%zl(x%Nlvgrid),x%zu(x%Nlvgrid))
    allocate(x%mz(x%Nlevels),x%mzl(x%Nlevels),x%mzu(x%Nlevels))
    
    ! --- Model vertical levels ---
    ! Use height levels of first model gridbox
    x%mz  = gbx%zlev(1,:)
    x%mzl = gbx%zlev_half(1,:)
    x%mzu(1:x%Nlevels-1) = gbx%zlev_half(1,2:x%Nlevels)
    x%mzu(x%Nlevels) = gbx%zlev(1,x%Nlevels) + (gbx%zlev(1,x%Nlevels) - x%mzl(x%Nlevels))
    
    if (use_vgrid) then
      ! --- Initialise to zero ---
      x%z  = 0.0
      x%zl = 0.0
      x%zu = 0.0
      if (cloudsat) then ! --- CloudSat grid requested ---
         zstep = 480.0
      else
         ! Other grid requested. Constant vertical spacing with top at 20 km
         zstep = 20000.0/x%Nlvgrid
      endif
      do i=1,x%Nlvgrid
         x%zl(i) = (i-1)*zstep
         x%zu(i) = i*zstep
      enddo
      x%z = (x%zl + x%zu)/2.0
    else
      x%z  = x%mz
      x%zl = x%mzl
      x%zu = x%mzu
    endif
    
  END SUBROUTINE CONSTRUCT_COSP_VGRID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------ SUBROUTINE FREE_COSP_VGRID ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_VGRID(x)
    type(cosp_vgrid),intent(inout) :: x

    deallocate(x%z, x%zl, x%zu, x%mz, x%mzl, x%mzu)
  END SUBROUTINE FREE_COSP_VGRID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_SGLIDAR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_SGLIDAR(cfg,Npoints,Ncolumns,Nlevels,Nhydro,Nrefl,x)
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Ncolumns ! Number of subgrid columns
    integer,intent(in) :: Nlevels  ! Number of model levels
    integer,intent(in) :: Nhydro   ! Number of hydrometeors
    integer,intent(in) :: Nrefl    ! Number of parasol reflectances ! parasol
    type(cosp_sglidar),intent(out) :: x
    ! Local variables
    integer :: i,j,k,l,m
    
    ! Allocate minumum storage if simulator not used
    if (cfg%Llidar_sim) then
      i = Npoints
      j = Ncolumns
      k = Nlevels
      l = Nhydro
      m = Nrefl
    else
      i = 1
      j = 1
      k = 1
      l = 1
      m = 1
    endif
    
    ! Dimensions
    x%Npoints  = i
    x%Ncolumns = j
    x%Nlevels  = k
    x%Nhydro   = l
    x%Nrefl    = m
    
    ! --- Allocate arrays ---
    allocate(x%beta_mol(i,k), x%beta_tot(i,j,k), &
             x%tau_tot(i,j,k),x%refl(i,j,m))
    ! --- Initialise to zero ---
    x%beta_mol   = 0.0
    x%beta_tot   = 0.0
    x%tau_tot    = 0.0
    x%refl       = 0.0 ! parasol
  END SUBROUTINE CONSTRUCT_COSP_SGLIDAR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------ SUBROUTINE FREE_COSP_SGLIDAR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_SGLIDAR(x)
    type(cosp_sglidar),intent(inout) :: x

    deallocate(x%beta_mol, x%beta_tot, x%tau_tot, x%refl)
  END SUBROUTINE FREE_COSP_SGLIDAR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_SGRADAR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_SGRADAR(cfg,Npoints,Ncolumns,Nlevels,Nhydro,x)
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Ncolumns ! Number of subgrid columns
    integer,intent(in) :: Nlevels  ! Number of model levels
    integer,intent(in) :: Nhydro   ! Number of hydrometeors
    type(cosp_sgradar),intent(out) :: x
    ! Local variables
    integer :: i,j,k,l
    
    if (cfg%Lradar_sim) then
      i = Npoints
      j = Ncolumns
      k = Nlevels
      l = Nhydro
    else ! Allocate minumum storage if simulator not used
      i = 1
      j = 1
      k = 1
      l = 1
    endif
    
    ! Dimensions
    x%Npoints  = i
    x%Ncolumns = j
    x%Nlevels  = k
    x%Nhydro   = l
    
    ! --- Allocate arrays ---
    allocate(x%att_gas(i,k), x%Ze_tot(i,j,k))
    ! --- Initialise to zero ---
    x%att_gas   = 0.0
    x%Ze_tot    = 0.0
    ! The following line give a compilation error on the Met Office NEC
!     call zero_real(x%Z_hydro, x%att_hydro)
!     f90: error(666): cosp_types.f90, line nnn:
!                                        Actual argument corresponding to dummy
!                                        argument of ELEMENTAL subroutine
!                                        "zero_real" with INTENET(OUT) attribute
!                                        is not array.
  END SUBROUTINE CONSTRUCT_COSP_SGRADAR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------ SUBROUTINE FREE_COSP_SGRADAR ----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_SGRADAR(x)
    type(cosp_sgradar),intent(inout) :: x

    deallocate(x%att_gas, x%Ze_tot)
  END SUBROUTINE FREE_COSP_SGRADAR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------- SUBROUTINE CONSTRUCT_COSP_RADARSTATS ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_RADARSTATS(cfg,Npoints,Ncolumns,Nlevels,Nhydro,x)
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Ncolumns ! Number of subgrid columns
    integer,intent(in) :: Nlevels  ! Number of model levels
    integer,intent(in) :: Nhydro   ! Number of hydrometeors
    type(cosp_radarstats),intent(out) :: x    
    ! Local variables
    integer :: i,j,k,l
    
    ! Allocate minumum storage if simulator not used
    if (cfg%Lradar_sim) then
      i = Npoints
      j = Ncolumns
      k = Nlevels
      l = Nhydro
    else
      i = 1
      j = 1
      k = 1
      l = 1
    endif
    
    ! Dimensions
    x%Npoints  = i
    x%Ncolumns = j
    x%Nlevels  = k
    x%Nhydro   = l
    
    ! --- Allocate arrays ---
    allocate(x%cfad_ze(i,DBZE_BINS,k),x%lidar_only_freq_cloud(i,k))
    allocate(x%radar_lidar_tcc(i),x%radar_tcc(i),x%radar_tcc_2(i)) !modified by ZYY
    ! --- Initialise to zero ---
    x%cfad_ze = 0.0
    x%lidar_only_freq_cloud = 0.0
    x%radar_lidar_tcc = 0.0
    x%radar_tcc = 0.0 !modified by ZYY
    x%radar_tcc_2 = 0.0 !modified by ZYY
  END SUBROUTINE CONSTRUCT_COSP_RADARSTATS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------ SUBROUTINE FREE_COSP_RADARSTATS -------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_RADARSTATS(x)
    type(cosp_radarstats),intent(inout) :: x

    deallocate(x%cfad_ze,x%lidar_only_freq_cloud,x%radar_lidar_tcc &
              ,x%radar_tcc,x%radar_tcc_2) !modified by ZYY
  END SUBROUTINE FREE_COSP_RADARSTATS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------- SUBROUTINE CONSTRUCT_COSP_LIDARSTATS ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_LIDARSTATS(cfg,Npoints,Ncolumns,Nlevels,Nhydro,Nrefl,x)
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Ncolumns ! Number of subgrid columns
    integer,intent(in) :: Nlevels  ! Number of model levels
    integer,intent(in) :: Nhydro   ! Number of hydrometeors
    integer,intent(in) :: Nrefl    ! Number of parasol reflectance
    type(cosp_lidarstats),intent(out) :: x
    ! Local variables
    integer :: i,j,k,l,m
    
    ! Allocate minumum storage if simulator not used
    if (cfg%Llidar_sim) then
      i = Npoints
      j = Ncolumns
      k = Nlevels
      l = Nhydro
      m = Nrefl
    else
      i = 1
      j = 1
      k = 1
      l = 1
      m = 1
    endif
    
    ! Dimensions
    x%Npoints  = i
    x%Ncolumns = j
    x%Nlevels  = k
    x%Nhydro   = l
    x%Nrefl    = m
    
    ! --- Allocate arrays ---
    allocate(x%srbval(SR_BINS),x%cfad_sr(i,SR_BINS,k), & 
             x%lidarcld(i,k), x%cldlayer(i,LIDAR_NCAT), x%parasolrefl(i,m))
    ! --- Initialise to zero ---
    x%srbval    = 0.0
    x%cfad_sr   = 0.0
    x%lidarcld  = 0.0
    x%cldlayer  = 0.0
    x%parasolrefl  = 0.0
  END SUBROUTINE CONSTRUCT_COSP_LIDARSTATS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------ SUBROUTINE FREE_COSP_LIDARSTATS -------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_LIDARSTATS(x)
    type(cosp_lidarstats),intent(inout) :: x

    deallocate(x%srbval, x%cfad_sr, x%lidarcld, x%cldlayer, x%parasolrefl)
  END SUBROUTINE FREE_COSP_LIDARSTATS
 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_SUBGRID ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_SUBGRID(Npoints,Ncolumns,Nlevels,y)
    integer,intent(in) :: Npoints, & ! Number of gridpoints
                                        Ncolumns, & ! Number of columns
                                        Nlevels   ! Number of levels
    type(cosp_subgrid),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels

    ! --- Allocate arrays ---
    allocate(y%frac_out(Npoints,Ncolumns,Nlevels))
    if (Ncolumns > 1) then
      allocate(y%prec_frac(Npoints,Ncolumns,Nlevels))
    else ! CRM mode, not needed
      allocate(y%prec_frac(1,1,1))
    endif
    ! --- Initialise to zero ---
    y%prec_frac = 0.0
    y%frac_out  = 0.0
    ! The following line gives a compilation error on the Met Office NEC
!     call zero_real(y%mr_hydro)
!     f90: error(666): cosp_types.f90, line nnn:
!                                        Actual argument corresponding to dummy
!                                        argument of ELEMENTAL subroutine
!                                        "zero_real" with INTENET(OUT) attribute
!                                        is not array.

  END SUBROUTINE CONSTRUCT_COSP_SUBGRID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_SUBGRID -----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_SUBGRID(y)
    type(cosp_subgrid),intent(inout) :: y
    
    ! --- Deallocate arrays ---
    deallocate(y%prec_frac, y%frac_out)
        
  END SUBROUTINE FREE_COSP_SUBGRID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_SGHYDRO -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_SGHYDRO(Npoints,Ncolumns,Nlevels,Nhydro,y)
    integer,intent(in) :: Npoints, & ! Number of gridpoints
                                        Ncolumns, & ! Number of columns
                                        Nhydro, & ! Number of hydrometeors
                                        Nlevels   ! Number of levels
    type(cosp_sghydro),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
    y%Nhydro   = Nhydro

    ! --- Allocate arrays ---
    allocate(y%mr_hydro(Npoints,Ncolumns,Nlevels,Nhydro), &
             y%Reff(Npoints,Ncolumns,Nlevels,Nhydro))
    ! --- Initialise to zero ---
    y%mr_hydro = 0.0
    y%Reff     = 0.0

  END SUBROUTINE CONSTRUCT_COSP_SGHYDRO

 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_SGHYDRO -----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_SGHYDRO(y)
    type(cosp_sghydro),intent(inout) :: y
    
    ! --- Deallocate arrays ---
    deallocate(y%mr_hydro, y%Reff)
        
  END SUBROUTINE FREE_COSP_SGHYDRO
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_GRIDBOX ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_GRIDBOX(time,time_bnds,radar_freq,surface_radar,use_mie_tables,use_gas_abs,do_ray,melt_lay,k2, &
                                   Npoints,Nlevels,Ncolumns,Nhydro,Nprmts_max_hydro,Naero,Nprmts_max_aero,Npoints_it, & 
                                   lidar_ice_type,isccp_top_height,isccp_top_height_direction,isccp_overlap,isccp_emsfc_lw, &
                                   use_precipitation_fluxes,use_reff, &
                                   ! RTTOV inputs
                                   Plat,Sat,Inst,Nchan,ZenAng,Ichan,SurfEm,co2,ch4,n2o,co,&
                                   y)
!jsb double precision -> real*8
    real*8,intent(in) :: time ! Time since start of run [days] 
    real*8,intent(in) :: time_bnds(2) ! Time boundaries
! jsb
    real,intent(in)    :: radar_freq, & ! Radar frequency [GHz]
                          k2            ! |K|^2, -1=use frequency dependent default
    integer,intent(in) :: &
        surface_radar, &  ! surface=1,spaceborne=0
        use_mie_tables, & ! use a precomputed lookup table? yes=1,no=0,2=use first column everywhere
        use_gas_abs, &    ! include gaseous absorption? yes=1,no=0
        do_ray, &         ! calculate/output Rayleigh refl=1, not=0
        melt_lay          ! melting layer model off=0, on=1
    integer,intent(in) :: Npoints   ! Number of gridpoints
    integer,intent(in) :: Nlevels   ! Number of levels
    integer,intent(in) :: Ncolumns  ! Number of columns
    integer,intent(in) :: Nhydro    ! Number of hydrometeors
    integer,intent(in) :: Nprmts_max_hydro    ! Max number of parameters for hydrometeor size distributions
    integer,intent(in) :: Naero    ! Number of aerosol species
    integer,intent(in) :: Nprmts_max_aero    ! Max number of parameters for aerosol size distributions
    integer,intent(in) :: Npoints_it   ! Number of gridpoints processed in one iteration
    integer,intent(in) :: lidar_ice_type ! Ice particle shape in lidar calculations (0=ice-spheres ; 1=ice-non-spherical)
    integer,intent(in) :: isccp_top_height
    integer,intent(in) :: isccp_top_height_direction
    integer,intent(in) :: isccp_overlap
    real,intent(in)    :: isccp_emsfc_lw
    logical,intent(in) :: use_precipitation_fluxes,use_reff
    integer,intent(in) :: Plat
    integer,intent(in) :: Sat
    integer,intent(in) :: Inst
    integer,intent(in) :: Nchan
    integer,intent(in) :: Ichan(Nchan)
    real,intent(in)    :: SurfEm(Nchan)
    real,intent(in)    :: ZenAng
    real,intent(in)    :: co2,ch4,n2o,co
    type(cosp_gridbox),intent(out) :: y

        
    ! local variables
    integer i, cnt_ice, cnt_liq
    real*8  :: delt, deltp
 
    ! Dimensions and scalars
    y%radar_freq       = radar_freq
    y%surface_radar    = surface_radar
    y%use_mie_tables   = use_mie_tables
    y%use_gas_abs      = use_gas_abs
    y%do_ray           = do_ray
    y%melt_lay         = melt_lay
    y%k2               = k2
    y%Npoints          = Npoints
    y%Nlevels          = Nlevels
    y%Ncolumns         = Ncolumns
    y%Nhydro           = Nhydro
    y%Nprmts_max_hydro = Nprmts_max_hydro
    y%Naero            = Naero
    y%Nprmts_max_aero  = Nprmts_max_aero
    y%Npoints_it       = Npoints_it
    y%lidar_ice_type   = lidar_ice_type
    y%isccp_top_height = isccp_top_height
    y%isccp_top_height_direction = isccp_top_height_direction
    y%isccp_overlap    = isccp_overlap
    y%isccp_emsfc_lw   = isccp_emsfc_lw
    y%use_precipitation_fluxes = use_precipitation_fluxes
    y%use_reff = use_reff
    
    y%time      = time
    y%time_bnds = time_bnds
    
    ! RTTOV parameters
    y%Plat   = Plat
    y%Sat    = Sat
    y%Inst   = Inst
    y%Nchan  = Nchan
    y%ZenAng = ZenAng
    y%co2    = co2
    y%ch4    = ch4
    y%n2o    = n2o
    y%co     = co

    ! --- Allocate arrays ---
    ! Gridbox information (Npoints,Nlevels)
    allocate(y%zlev(Npoints,Nlevels), y%zlev_half(Npoints,Nlevels), y%dlev(Npoints,Nlevels), &
             y%p(Npoints,Nlevels), y%ph(Npoints,Nlevels), y%T(Npoints,Nlevels), &
             y%q(Npoints,Nlevels), y%sh(Npoints,Nlevels), &
             !jsb
             y%dtau_s(Npoints,Nlevels), y%dtau_c(Npoints,Nlevels),y%dtau_s_snow(Npoints,Nlevels), &
             y%dem_s(Npoints,Nlevels), y%dem_c(Npoints,Nlevels),y%dem_s_snow(Npoints,Nlevels), &
             !jsb
             y%tca(Npoints,Nlevels), y%cca(Npoints,Nlevels), &
             y%rain_ls(Npoints,Nlevels), y%rain_cv(Npoints,Nlevels), y%grpl_ls(Npoints,Nlevels), &
             y%snow_ls(Npoints,Nlevels), y%snow_cv(Npoints,Nlevels),y%mr_ozone(Npoints,Nlevels))
             
             
    ! Surface information and geolocation (Npoints)
    allocate(y%longitude(Npoints),y%latitude(Npoints),y%psfc(Npoints), y%land(Npoints), &
             y%sunlit(Npoints),y%skt(Npoints),y%sfc_height(Npoints),y%u_wind(Npoints),y%v_wind(Npoints))
    ! Hydrometeors concentration and distribution parameters
    allocate(y%mr_hydro(Npoints,Nlevels,Nhydro), &
             y%dist_prmts_hydro(Nprmts_max_hydro,Nhydro), &
             y%Reff(Npoints,Nlevels,Nhydro))
    ! Aerosols concentration and distribution parameters
    allocate(y%conc_aero(Npoints,Nlevels,Naero), y%dist_type_aero(Naero), &
             y%dist_prmts_aero(Npoints,Nlevels,Nprmts_max_aero,Naero))
    
    ! RTTOV channels and sfc. emissivity
    allocate(y%ichan(Nchan),y%surfem(Nchan))
    
    ! RTTOV parameters
    y%ichan   =  ichan
    y%surfem  =  surfem
    
    ! --- Initialise to zero ---
    y%zlev      = 0.0
    y%zlev_half = 0.0
    y%dlev      = 0.0
    y%p         = 0.0
    y%ph        = 0.0
    y%T         = 0.0
    y%q         = 0.0
    y%sh        = 0.0
    y%dtau_s    = 0.0
    y%dtau_c    = 0.0
    ! jsb
    y%dtau_s_snow = 0.0
    ! jsb
    y%dem_s     = 0.0
    y%dem_c     = 0.0
    ! jsb
    y%dem_s_snow = 0.0
    ! jsb
    y%tca       = 0.0
    y%cca       = 0.0
    y%rain_ls   = 0.0
    y%rain_cv   = 0.0
    y%grpl_ls   = 0.0
    y%snow_ls   = 0.0
    y%snow_cv   = 0.0
    y%Reff      = 0.0
    y%mr_ozone  = 0.0
    y%u_wind    = 0.0
    y%v_wind    = 0.0

    
    ! (Npoints)
!     call zero_real(y%psfc, y%land)
    y%longitude = 0.0
    y%latitude = 0.0
    y%psfc = 0.0
    y%land = 0.0
    y%sunlit = 0.0
    y%skt = 0.0
    y%sfc_height = 0.0
    ! (Npoints,Nlevels,Nhydro)
!     y%fr_hydro = 0.0
    y%mr_hydro = 0.0
    ! Others
    y%dist_prmts_hydro = 0.0 ! (Nprmts_max_hydro,Nhydro)
    y%conc_aero        = 0.0 ! (Npoints,Nlevels,Naero)
    y%dist_type_aero   = 0   ! (Naero)
    y%dist_prmts_aero  = 0.0 ! (Npoints,Nlevels,Nprmts_max_aero,Naero)

    y%hp%p1 = 0.0
    y%hp%p2 = 0.0
    y%hp%p3 = 0.0
    y%hp%dmin = 0.0
    y%hp%dmax = 0.0
    y%hp%apm = 0.0
    y%hp%bpm = 0.0
    y%hp%rho = 0.0
    y%hp%dtype = 0
    y%hp%col = 0
    y%hp%cp = 0
    y%hp%phase = 0
    y%hp%scaled = .false.
    y%hp%z_flag = .false.
    y%hp%Ze_scaled = 0.0
    y%hp%Zr_scaled = 0.0
    y%hp%kr_scaled = 0.0
    y%hp%fc = 0.0
    y%hp%rho_eff = 0.0
    y%hp%ifc = 0
    y%hp%idd = 0
    y%mt%freq = 0.0
    y%mt%tt = 0.0
    y%mt%f = 0.0
    y%mt%D = 0.0
    y%mt%qext = 0.0
    y%mt%qbsca = 0.0
    y%mt%phase = 0
    
    
    ! --- Initialize the distributional parameters for hydrometeors
    y%dist_prmts_hydro( 1,:) = HCLASS_TYPE(:)
    y%dist_prmts_hydro( 2,:) = HCLASS_COL(:)
    y%dist_prmts_hydro( 3,:) = HCLASS_PHASE(:)
    y%dist_prmts_hydro( 4,:) = HCLASS_CP(:)
    y%dist_prmts_hydro( 5,:) = HCLASS_DMIN(:)
    y%dist_prmts_hydro( 6,:) = HCLASS_DMAX(:)
    y%dist_prmts_hydro( 7,:) = HCLASS_APM(:)
    y%dist_prmts_hydro( 8,:) = HCLASS_BPM(:)
    y%dist_prmts_hydro( 9,:) = HCLASS_RHO(:)
    y%dist_prmts_hydro(10,:) = HCLASS_P1(:)
    y%dist_prmts_hydro(11,:) = HCLASS_P2(:)
    y%dist_prmts_hydro(12,:) = HCLASS_P3(:)

    ! the following code added by roj to initialize structures used by radar simulator, Feb 2008
    call load_hydrometeor_classes(y%Nprmts_max_hydro,y%dist_prmts_hydro(:,:),y%hp,y%Nhydro)

    ! load mie tables ?
    if (y%use_mie_tables == 1) then
      print *, '%%% COSP: Mie tables option for Quickbem not supported'
      stop
!         ! ----- Mie tables ----
!           mie_table_name='mie_table.dat'
!         call load_mie_table(mie_table_name,y%mt)
!   
!       !   :: D specified by table ... not must match that used when mie LUT generated!
!       y%nsizes = mt_nd
!       allocate(y%D(y%nsizes))
!       y%D = y%mt%D

    else
       ! otherwise we still need to initialize temperature arrays for Ze scaling (which is only done when not using mie table)
       
       cnt_ice=19
       cnt_liq=20
!        if (.not.(allocated(mt_ttl).and.allocated(mt_tti))) then
!           allocate(mt_ttl(cnt_liq),mt_tti(cnt_ice))  ! note needed as this is global array ... 
!                                                      ! which should be changed in the future
!        endif
          
       do i=1,cnt_ice
          mt_tti(i)=(i-1)*5-90
       enddo
    
       do i=1,cnt_liq
          mt_ttl(i)=(i-1)*5 - 60
       enddo 
    
       allocate(y%mt_ttl(cnt_liq),y%mt_tti(cnt_ice))

       y%mt_ttl = mt_ttl
       y%mt_tti = mt_tti

! !------ OLD code in v0.1 ---------------------------
!        allocate(mt_ttl(2),mt_tti(2))
!        allocate(y%mt_ttl(2),y%mt_tti(2))
!        mt_ttl = 0.0
!        mt_tti = 0.0
!        y%mt_ttl = mt_ttl
!        y%mt_tti = mt_tti
! !---------------------------------------------------
       
       ! :: D created on a log-linear scale
       y%nsizes = nd
       delt = (log(dmax)-log(dmin))/(y%nsizes-1)
       deltp = exp(delt)
       allocate(y%D(y%nsizes))
       y%D(1) = dmin
       do i=2,y%nsizes
          y%D(i) = y%D(i-1)*deltp
       enddo   
   
    endif


END SUBROUTINE CONSTRUCT_COSP_GRIDBOX

  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_GRIDBOX -----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_GRIDBOX(y,dglobal)
    type(cosp_gridbox),intent(inout) :: y
    logical,intent(in),optional :: dglobal

    ! --- Free arrays ---
    deallocate(y%D,y%mt_ttl,y%mt_tti)   ! added by roj Feb 2008
!     if (.not.present(dglobal)) deallocate(mt_ttl,mt_tti)
    
!     deallocate(y%hp%p1,y%hp%p2,y%hp%p3,y%hp%dmin,y%hp%dmax,y%hp%apm,y%hp%bpm,y%hp%rho, &
!               y%hp%dtype,y%hp%col,y%hp%cp,y%hp%phase,y%hp%scaled, &
!               y%hp%z_flag,y%hp%Ze_scaled,y%hp%Zr_scaled,y%hp%kr_scaled, &
!               y%hp%fc, y%hp%rho_eff, y%hp%ifc, y%hp%idd)
!     deallocate(y%mt%freq, y%mt%tt, y%mt%f, y%mt%D, y%mt%qext, y%mt%qbsca, y%mt%phase)
    
    deallocate(y%zlev, y%zlev_half, y%dlev, y%p, y%ph, y%T, y%q, &
               y%sh, y%dtau_s, y%dtau_c, y%dem_s, y%dem_c, &
               !jsb
               y%dtau_s_snow,y%dem_s_snow, &
               !jsb
               y%longitude,y%latitude,y%psfc, y%land, y%tca, y%cca, &
               y%mr_hydro, y%dist_prmts_hydro, &
               y%conc_aero, y%dist_type_aero, y%dist_prmts_aero, &
               y%rain_ls, y%rain_cv, y%snow_ls, y%snow_cv, y%grpl_ls, &
               y%sunlit, y%skt, y%sfc_height, y%Reff,y%ichan,y%surfem, &
               y%mr_ozone,y%u_wind,y%v_wind)
 
  END SUBROUTINE FREE_COSP_GRIDBOX
  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_GRIDBOX_CPHP ----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_GRIDBOX_CPHP(x,y)
    type(cosp_gridbox),intent(in) :: x
    type(cosp_gridbox),intent(inout) :: y
    
    integer :: i,j,k,sz(3)
    ! jsb
    real*8 :: tny
    ! jsb
    
    tny = tiny(tny)
    y%hp%p1      = x%hp%p1
    y%hp%p2      = x%hp%p2
    y%hp%p3      = x%hp%p3
    y%hp%dmin    = x%hp%dmin
    y%hp%dmax    = x%hp%dmax
    y%hp%apm     = x%hp%apm
    y%hp%bpm     = x%hp%bpm
    y%hp%rho     = x%hp%rho
    y%hp%dtype   = x%hp%dtype
    y%hp%col     = x%hp%col
    y%hp%cp      = x%hp%cp
    y%hp%phase   = x%hp%phase

    y%hp%fc      = x%hp%fc
    y%hp%rho_eff = x%hp%rho_eff
    y%hp%ifc     = x%hp%ifc
    y%hp%idd     = x%hp%idd
    sz = shape(x%hp%z_flag)
    do k=1,sz(3)
      do j=1,sz(2)
        do i=1,sz(1)
           if (x%hp%scaled(i,k))   y%hp%scaled(i,k)      = .true.
           if (x%hp%z_flag(i,j,k)) y%hp%z_flag(i,j,k)    = .true.
           if (abs(x%hp%Ze_scaled(i,j,k)) > tny) y%hp%Ze_scaled(i,j,k) = x%hp%Ze_scaled(i,j,k)
           if (abs(x%hp%Zr_scaled(i,j,k)) > tny) y%hp%Zr_scaled(i,j,k) = x%hp%Zr_scaled(i,j,k)
           if (abs(x%hp%kr_scaled(i,j,k)) > tny) y%hp%kr_scaled(i,j,k) = x%hp%kr_scaled(i,j,k)
        enddo
      enddo
    enddo
    
END SUBROUTINE COSP_GRIDBOX_CPHP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_GRIDBOX_CPSECTION -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_GRIDBOX_CPSECTION(ix,iy,x,y)
    integer,intent(in),dimension(2) :: ix,iy
    type(cosp_gridbox),intent(in) :: x
    type(cosp_gridbox),intent(inout) :: y
    
    ! --- Copy arrays without Npoints as dimension ---
    y%dist_prmts_hydro = x%dist_prmts_hydro
    y%dist_type_aero   = x%dist_type_aero
    y%D                = x%D
    y%mt_ttl           = x%mt_ttl
    y%mt_tti           = x%mt_tti
    
    
!     call cosp_gridbox_cphp(x,y)    
    
    ! 1D
    y%longitude(iy(1):iy(2))  = x%longitude(ix(1):ix(2))
    y%latitude(iy(1):iy(2))   = x%latitude(ix(1):ix(2))
    y%psfc(iy(1):iy(2))       = x%psfc(ix(1):ix(2))
    y%land(iy(1):iy(2))       = x%land(ix(1):ix(2))
    y%sunlit(iy(1):iy(2))     = x%sunlit(ix(1):ix(2))
    y%skt(iy(1):iy(2))        = x%skt(ix(1):ix(2))
    y%sfc_height(iy(1):iy(2)) = x%sfc_height(ix(1):ix(2))
    y%u_wind(iy(1):iy(2))     = x%u_wind(ix(1):ix(2))
    y%v_wind(iy(1):iy(2))     = x%v_wind(ix(1):ix(2))
    ! 2D
    y%zlev(iy(1):iy(2),:)      = x%zlev(ix(1):ix(2),:)
    y%zlev_half(iy(1):iy(2),:) = x%zlev_half(ix(1):ix(2),:)
    y%dlev(iy(1):iy(2),:)      = x%dlev(ix(1):ix(2),:)
    y%p(iy(1):iy(2),:)         = x%p(ix(1):ix(2),:)
    y%ph(iy(1):iy(2),:)        = x%ph(ix(1):ix(2),:)
    y%T(iy(1):iy(2),:)         = x%T(ix(1):ix(2),:)
    y%q(iy(1):iy(2),:)         = x%q(ix(1):ix(2),:)
    y%sh(iy(1):iy(2),:)        = x%sh(ix(1):ix(2),:)
    y%dtau_s(iy(1):iy(2),:)    = x%dtau_s(ix(1):ix(2),:)
    y%dtau_c(iy(1):iy(2),:)    = x%dtau_c(ix(1):ix(2),:)
    !jsb
    y%dtau_s_snow(iy(1):iy(2),:) = x%dtau_s_snow(ix(1):ix(2),:)
    !jsb
    y%dem_s(iy(1):iy(2),:)     = x%dem_s(ix(1):ix(2),:)
    y%dem_c(iy(1):iy(2),:)     = x%dem_c(ix(1):ix(2),:)
    !jsb
    y%dem_s_snow(iy(1):iy(2),:) = x%dem_s_snow(ix(1):ix(2),:)
    !jsb
    y%tca(iy(1):iy(2),:)       = x%tca(ix(1):ix(2),:)
    y%cca(iy(1):iy(2),:)       = x%cca(ix(1):ix(2),:)
    y%rain_ls(iy(1):iy(2),:)   = x%rain_ls(ix(1):ix(2),:)
    y%rain_cv(iy(1):iy(2),:)   = x%rain_cv(ix(1):ix(2),:)
    y%grpl_ls(iy(1):iy(2),:)   = x%grpl_ls(ix(1):ix(2),:)
    y%snow_ls(iy(1):iy(2),:)   = x%snow_ls(ix(1):ix(2),:)
    y%snow_cv(iy(1):iy(2),:)   = x%snow_cv(ix(1):ix(2),:)
    y%mr_ozone(iy(1):iy(2),:)  = x%mr_ozone(ix(1):ix(2),:)
    ! 3D
    y%Reff(iy(1):iy(2),:,:)      = x%Reff(ix(1):ix(2),:,:)
    y%conc_aero(iy(1):iy(2),:,:) = x%conc_aero(ix(1):ix(2),:,:)
    y%mr_hydro(iy(1):iy(2),:,:)  = x%mr_hydro(ix(1):ix(2),:,:)
    ! 4D
    y%dist_prmts_aero(iy(1):iy(2),:,:,:) = x%dist_prmts_aero(ix(1):ix(2),:,:,:)

END SUBROUTINE COSP_GRIDBOX_CPSECTION
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_SUBGRID_CPSECTION -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_SUBGRID_CPSECTION(ix,iy,x,y)
    integer,intent(in),dimension(2) :: ix,iy
    type(cosp_subgrid),intent(in) :: x
    type(cosp_subgrid),intent(inout) :: y
    
    y%prec_frac(iy(1):iy(2),:,:)  = x%prec_frac(ix(1):ix(2),:,:)
    y%frac_out(iy(1):iy(2),:,:)   = x%frac_out(ix(1):ix(2),:,:)
END SUBROUTINE COSP_SUBGRID_CPSECTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_SGRADAR_CPSECTION -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_SGRADAR_CPSECTION(ix,iy,x,y)
    integer,intent(in),dimension(2) :: ix,iy
    type(cosp_sgradar),intent(in) :: x
    type(cosp_sgradar),intent(inout) :: y
    
    y%att_gas(iy(1):iy(2),:)  = x%att_gas(ix(1):ix(2),:)
    y%Ze_tot(iy(1):iy(2),:,:) = x%Ze_tot(ix(1):ix(2),:,:)
END SUBROUTINE COSP_SGRADAR_CPSECTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_SGLIDAR_CPSECTION -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_SGLIDAR_CPSECTION(ix,iy,x,y)
    integer,intent(in),dimension(2) :: ix,iy
    type(cosp_sglidar),intent(in) :: x
    type(cosp_sglidar),intent(inout) :: y
    
    y%beta_mol(iy(1):iy(2),:)       = x%beta_mol(ix(1):ix(2),:)
    y%beta_tot(iy(1):iy(2),:,:)     = x%beta_tot(ix(1):ix(2),:,:)
    y%tau_tot(iy(1):iy(2),:,:)      = x%tau_tot(ix(1):ix(2),:,:)
    y%refl(iy(1):iy(2),:,:)         = x%refl(ix(1):ix(2),:,:)
END SUBROUTINE COSP_SGLIDAR_CPSECTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_ISCCP_CPSECTION -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_ISCCP_CPSECTION(ix,iy,x,y)
    integer,intent(in),dimension(2) :: ix,iy
    type(cosp_isccp),intent(in) :: x
    type(cosp_isccp),intent(inout) :: y
            
    y%fq_isccp(iy(1):iy(2),:,:)  = x%fq_isccp(ix(1):ix(2),:,:)
    y%totalcldarea(iy(1):iy(2))  = x%totalcldarea(ix(1):ix(2))
    y%meantb(iy(1):iy(2))        = x%meantb(ix(1):ix(2))
    y%meantbclr(iy(1):iy(2))     = x%meantbclr(ix(1):ix(2))
    y%meanptop(iy(1):iy(2))      = x%meanptop(ix(1):ix(2))
    y%meantaucld(iy(1):iy(2))    = x%meantaucld(ix(1):ix(2))
    y%meanalbedocld(iy(1):iy(2)) = x%meanalbedocld(ix(1):ix(2))
    y%boxtau(iy(1):iy(2),:)      = x%boxtau(ix(1):ix(2),:)
    y%boxptop(iy(1):iy(2),:)     = x%boxptop(ix(1):ix(2),:)
END SUBROUTINE COSP_ISCCP_CPSECTION


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_MISR_CPSECTION -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_MISR_CPSECTION(ix,iy,x,y)
    integer,intent(in),dimension(2) :: ix,iy
    type(cosp_misr),intent(in) :: x
    type(cosp_misr),intent(inout) :: y
            
    y%fq_MISR(iy(1):iy(2),:,:)                 = x%fq_MISR(ix(1):ix(2),:,:)
    y%MISR_meanztop(iy(1):iy(2))               = x%MISR_meanztop(ix(1):ix(2))
    y%MISR_cldarea(iy(1):iy(2))                = x%MISR_cldarea(ix(1):ix(2))
    y%MISR_dist_model_layertops(iy(1):iy(2),:) = x%MISR_dist_model_layertops(ix(1):ix(2),:)
END SUBROUTINE COSP_MISR_CPSECTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_RTTOV_CPSECTION -------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_RTTOV_CPSECTION(ix,iy,x,y)
    integer,intent(in),dimension(2) :: ix,iy
    type(cosp_rttov),intent(in) :: x
    type(cosp_rttov),intent(inout) :: y
            
    y%tbs(iy(1):iy(2),:) = x%tbs(ix(1):ix(2),:)
END SUBROUTINE COSP_RTTOV_CPSECTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_RADARSTATS_CPSECTION --------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_RADARSTATS_CPSECTION(ix,iy,x,y)
    integer,intent(in),dimension(2) :: ix,iy
    type(cosp_radarstats),intent(in) :: x
    type(cosp_radarstats),intent(inout) :: y
            
    y%cfad_ze(iy(1):iy(2),:,:)             = x%cfad_ze(ix(1):ix(2),:,:)
    y%radar_lidar_tcc(iy(1):iy(2))         = x%radar_lidar_tcc(ix(1):ix(2))
    y%radar_tcc(iy(1):iy(2))         = x%radar_tcc(ix(1):ix(2)) !modified by ZYY
    y%radar_tcc_2(iy(1):iy(2))         = x%radar_tcc_2(ix(1):ix(2)) !modified by ZYY
    y%lidar_only_freq_cloud(iy(1):iy(2),:) = x%lidar_only_freq_cloud(ix(1):ix(2),:)
END SUBROUTINE COSP_RADARSTATS_CPSECTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_LIDARSTATS_CPSECTION --------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_LIDARSTATS_CPSECTION(ix,iy,x,y)
    integer,intent(in),dimension(2) :: ix,iy
    type(cosp_lidarstats),intent(in) :: x
    type(cosp_lidarstats),intent(inout) :: y
            
    y%srbval                     = x%srbval
    y%cfad_sr(iy(1):iy(2),:,:)   = x%cfad_sr(ix(1):ix(2),:,:)
    y%lidarcld(iy(1):iy(2),:)    = x%lidarcld(ix(1):ix(2),:)
    y%cldlayer(iy(1):iy(2),:)    = x%cldlayer(ix(1):ix(2),:)
    y%parasolrefl(iy(1):iy(2),:) = x%parasolrefl(ix(1):ix(2),:)
END SUBROUTINE COSP_LIDARSTATS_CPSECTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- PRINT SUBROUTINES --------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_GRIDBOX_PRINT(x)
    type(cosp_gridbox),intent(in) :: x

    print *, '%%%%----- Information on COSP_GRIDBOX ------'
    ! Scalars and dimensions
    print *,  x%Npoints
    print *,  x%Nlevels
    print *,  x%Ncolumns
    print *,  x%Nhydro
    print *,  x%Nprmts_max_hydro
    print *,  x%Naero
    print *,  x%Nprmts_max_aero
    print *,  x%Npoints_it
    
    ! Time [days]
    print *,  x%time
    
    ! Radar ancillary info
    print *,  x%radar_freq, &
            x%k2
    print *,  x%surface_radar, &
           x%use_mie_tables, &
           x%use_gas_abs, &
           x%do_ray, &
           x%melt_lay

!               print *,  'shape(x%): ',shape(x%)
 
!     type(class_param) ::  hp  ! structure used by radar simulator to store Ze and N scaling constants and other information
!     type(mie)::  mt       ! structure used by radar simulator to store mie LUT information
    print *,  x%nsizes
    print *,  'shape(x%D): ',shape(x%D)
    print *,  'shape(x%mt_ttl): ',shape(x%mt_ttl)
    print *,  'shape(x%mt_tti): ',shape(x%mt_tti)
    
    ! Lidar
    print *,  x%lidar_ice_type
    
    ! Radar
    print *,  x%use_precipitation_fluxes
    print *,  x%use_reff
    
    ! Geolocation (Npoints)
    print *,  'shape(x%longitude): ',shape(x%longitude)
    print *,  'shape(x%latitude): ',shape(x%latitude)
    ! Gridbox information (Npoints,Nlevels)
    print *,  'shape(x%zlev): ',shape(x%zlev)
    print *,  'shape(x%zlev_half): ',shape(x%zlev_half)
    print *,  'shape(x%dlev): ',shape(x%dlev)
    print *,  'shape(x%p): ',shape(x%p)
    print *,  'shape(x%ph): ',shape(x%ph)
    print *,  'shape(x%T): ',shape(x%T)
    print *,  'shape(x%q): ',shape(x%q)
    print *,  'shape(x%sh): ',shape(x%sh)
    print *,  'shape(x%dtau_s): ',shape(x%dtau_s)
    print *,  'shape(x%dtau_c): ',shape(x%dtau_c)
    print *,  'shape(x%dem_s): ',shape(x%dem_s)
    ! jsb
    print *,  'shape(x%dtau_s_snow): ',shape(x%dtau_s_snow)
    !jsb
    print *,  'shape(x%dem_c): ',shape(x%dem_c)
    ! jsb
    print *,  'shape(x%dem_s_snow): ',shape(x%dem_s_snow)
    !jsb
    print *,  'shape(x%mr_ozone): ',shape(x%mr_ozone)

    ! Point information (Npoints)
    print *,  'shape(x%land): ',shape(x%land)
    print *,  'shape(x%psfc): ',shape(x%psfc)
    print *,  'shape(x%sunlit): ',shape(x%sunlit)
    print *,  'shape(x%skt): ',shape(x%skt)
    print *,  'shape(x%sfc_height): ',shape(x%sfc_height)
    print *,  'shape(x%u_wind): ',shape(x%u_wind)
    print *,  'shape(x%v_wind): ',shape(x%v_wind)

    ! TOTAL and CONV cloud fraction for SCOPS
    print *,  'shape(x%tca): ',shape(x%tca)
    print *,  'shape(x%cca): ',shape(x%cca)
    ! Precipitation fluxes on model levels
    print *,  'shape(x%rain_ls): ',shape(x%rain_ls)
    print *,  'shape(x%rain_cv): ',shape(x%rain_cv)
    print *,  'shape(x%snow_ls): ',shape(x%snow_ls)
    print *,  'shape(x%snow_cv): ',shape(x%snow_cv)
    print *,  'shape(x%grpl_ls): ',shape(x%grpl_ls)
    ! Hydrometeors concentration and distribution parameters
    print *,  'shape(x%mr_hydro): ',shape(x%mr_hydro)
    print *,  'shape(x%dist_prmts_hydro): ',shape(x%dist_prmts_hydro)
    ! Effective radius [m]. (Npoints,Nlevels,Nhydro)
    print *,  'shape(x%Reff): ',shape(x%Reff)
    ! Aerosols concentration and distribution parameters
    print *,  'shape(x%conc_aero): ',shape(x%conc_aero)
    print *,  'shape(x%dist_type_aero): ',shape(x%dist_type_aero)
    print *,  'shape(x%dist_prmts_aero): ',shape(x%dist_prmts_aero)
    ! ISCCP simulator inputs
    print *, x%isccp_top_height
    print *, x%isccp_top_height_direction
    print *, x%isccp_overlap
    print *, x%isccp_emsfc_lw
  
    ! RTTOV inputs/options
    print *, x%plat
    print *, x%sat
    print *, x%inst
    print *, x%Nchan
    print *,  'shape(x%Ichan): ',x%Ichan
    print *,  'shape(x%Surfem): ',x%Surfem
    print *, x%ZenAng
    print *, x%co2,x%ch4,x%n2o,x%co
                
END SUBROUTINE COSP_GRIDBOX_PRINT

SUBROUTINE COSP_MISR_PRINT(x)
    type(cosp_misr),intent(in) :: x

    print *, '%%%%----- Information on COSP_MISR ------'
                
     ! Dimensions
    print *, x%Npoints
    print *, x%Ntau
    print *, x%Nlevels

     ! --- (npoints,ntau,nlevels)
     !  the fraction of the model grid box covered by each of the MISR cloud types
     print *,  'shape(x%fq_MISR): ',shape(x%fq_MISR)
     
     ! --- (npoints)
     print *,  'shape(x%MISR_meanztop): ',shape(x%MISR_meanztop)
     print *,  'shape(x%MISR_cldarea): ',shape(x%MISR_cldarea)
     ! --- (npoints,nlevels)
     print *,  'shape(x%MISR_dist_model_layertops): ',shape(x%MISR_dist_model_layertops)
    
END SUBROUTINE COSP_MISR_PRINT

SUBROUTINE COSP_ISCCP_PRINT(x)
    type(cosp_isccp),intent(in) :: x
            
    print *, x%Npoints
    print *, x%Ncolumns
    print *, x%Nlevels

    print *, '%%%%----- Information on COSP_ISCCP ------'
    
     print *, 'shape(x%fq_isccp): ',shape(x%fq_isccp)
     print *, 'shape(x%totalcldarea): ',shape(x%totalcldarea)
     print *, 'shape(x%meantb): ',shape(x%meantb)
     print *, 'shape(x%meantbclr): ',shape(x%meantbclr)
     
     print *, 'shape(x%meanptop): ',shape(x%meanptop)
     print *, 'shape(x%meantaucld): ',shape(x%meantaucld)
     print *, 'shape(x%meanalbedocld): ',shape(x%meanalbedocld)
     print *, 'shape(x%boxtau): ',shape(x%boxtau)
     print *, 'shape(x%boxptop): ',shape(x%boxptop)
END SUBROUTINE COSP_ISCCP_PRINT

SUBROUTINE COSP_VGRID_PRINT(x)
    type(cosp_vgrid),intent(in) :: x
            
    print *, '%%%%----- Information on COSP_VGRID ------'
    print *, x%use_vgrid
    print *, x%csat_vgrid
    print *, x%Npoints
    print *, x%Ncolumns
    print *, x%Nlevels
    print *, x%Nlvgrid
    ! Array with dimensions (Nlvgrid)
    print *, 'shape(x%z): ',shape(x%z)
    print *, 'shape(x%zl): ',shape(x%zl)
    print *, 'shape(x%zu): ',shape(x%zu)
    ! Array with dimensions (Nlevels)
    print *, 'shape(x%mz): ',shape(x%mz)
    print *, 'shape(x%mzl): ',shape(x%mzl)
    print *, 'shape(x%mzu): ',shape(x%mzu)
END SUBROUTINE COSP_VGRID_PRINT

SUBROUTINE COSP_SGLIDAR_PRINT(x)
    type(cosp_sglidar),intent(in) :: x
            
    print *, '%%%%----- Information on COSP_SGLIDAR ------'
    ! Dimensions
    print *, x%Npoints
    print *, x%Ncolumns
    print *, x%Nlevels
    print *, x%Nhydro
    print *, x%Nrefl
    ! Arrays with dimensions (Npoints,Nlevels)
    print *, 'shape(x%beta_mol): ',shape(x%beta_mol)
    ! Arrays with dimensions (Npoints,Ncolumns,Nlevels)
    print *, 'shape(x%beta_tot): ',shape(x%beta_tot)
    print *, 'shape(x%tau_tot): ',shape(x%tau_tot)
    ! Arrays with dimensions (Npoints,Ncolumns,Nrefl)
    print *, 'shape(x%refl): ',shape(x%refl)
END SUBROUTINE COSP_SGLIDAR_PRINT

SUBROUTINE COSP_SGRADAR_PRINT(x)
    type(cosp_sgradar),intent(in) :: x
            
    print *, '%%%%----- Information on COSP_SGRADAR ------'
    print *, x%Npoints
    print *, x%Ncolumns
    print *, x%Nlevels
    print *, x%Nhydro
    ! output vertical levels: spaceborne radar -> from TOA to SURFACE
    ! Arrays with dimensions (Npoints,Nlevels)
    print *, 'shape(x%att_gas): ', shape(x%att_gas)
    ! Arrays with dimensions (Npoints,Ncolumns,Nlevels)
    print *, 'shape(x%Ze_tot): ', shape(x%Ze_tot)
END SUBROUTINE COSP_SGRADAR_PRINT

SUBROUTINE COSP_RADARSTATS_PRINT(x)
    type(cosp_radarstats),intent(in) :: x
            
    print *, '%%%%----- Information on COSP_SGRADAR ------'
    print *, x%Npoints
    print *, x%Ncolumns
    print *, x%Nlevels
    print *, x%Nhydro
    print *, 'shape(x%cfad_ze): ',shape(x%cfad_ze)
    print *, 'shape(x%radar_lidar_tcc): ',shape(x%radar_lidar_tcc)
    print *, 'shape(x%radar_tcc): ',shape(x%radar_tcc) !modified by ZYY
    print *, 'shape(x%radar_tcc_2): ',shape(x%radar_tcc_2) !modified by ZYY
    print *, 'shape(x%lidar_only_freq_cloud): ',shape(x%lidar_only_freq_cloud)
END SUBROUTINE COSP_RADARSTATS_PRINT

SUBROUTINE COSP_LIDARSTATS_PRINT(x)
    type(cosp_lidarstats),intent(in) :: x
            
    print *, '%%%%----- Information on COSP_SGLIDAR ------'
    print *, x%Npoints
    print *, x%Ncolumns
    print *, x%Nlevels
    print *, x%Nhydro
    print *, x%Nrefl
    
    ! Arrays with dimensions (SR_BINS)
    print *, 'shape(x%srbval): ',shape(x%srbval)
    ! Arrays with dimensions (Npoints,SR_BINS,Nlevels)
    print *, 'shape(x%cfad_sr): ',shape(x%cfad_sr)
    ! Arrays with dimensions (Npoints,Nlevels)
    print *, 'shape(x%lidarcld): ',shape(x%lidarcld)
    ! Arrays with dimensions (Npoints,LIDAR_NCAT)
    print *, 'shape(x%cldlayer): ',shape(x%cldlayer)
    ! Arrays with dimensions (Npoints,PARASOL_NREFL)
    print *, 'shape(x%parasolrefl): ',shape(x%parasolrefl)
END SUBROUTINE COSP_LIDARSTATS_PRINT

SUBROUTINE COSP_SUBGRID_PRINT(x)
    type(cosp_subgrid),intent(in) :: x
            
    print *, '%%%%----- Information on COSP_SUBGRID ------'
    print *, x%Npoints
    print *, x%Ncolumns
    print *, x%Nlevels
    print *, x%Nhydro
    
    print *, 'shape(x%prec_frac): ',shape(x%prec_frac)
    print *, 'shape(x%frac_out): ',shape(x%frac_out)
END SUBROUTINE COSP_SUBGRID_PRINT

SUBROUTINE COSP_SGHYDRO_PRINT(x)
    type(cosp_sghydro),intent(in) :: x
            
    print *, '%%%%----- Information on COSP_SGHYDRO ------'
    print *, x%Npoints
    print *, x%Ncolumns
    print *, x%Nlevels
    print *, x%Nhydro
    
    print *, 'shape(x%mr_hydro): ',shape(x%mr_hydro)
    print *, 'shape(x%Reff): ',shape(x%Reff)
END SUBROUTINE COSP_SGHYDRO_PRINT

END MODULE MOD_COSP_TYPES
