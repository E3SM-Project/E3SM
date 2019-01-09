MODULE q3d_init_module
! Contains the programs related to initializing the channels

      USE shr_kind_mod,   only: r8 => shr_kind_r8
      USE vvm_data_types, only: channel_t

      USE parmsld, only: channel_l,channel_w,ntracer,nhalo,nhalo_adv, &
                         nk1,nk2,nk3,num_ext_layers,nk1_ext,nk2_ext
                         
      USE mphysics_variables, only: rho_int,rhol,dzl,pl0,pil0

      USE constld

      USE trace_gases, only: trace_gas_input

      USE rrtm_grid,   only: nrad_rad => nrad
      USE rrtm_vars,   only: pres,presi

      USE time_manager,         only: get_step_size     ! ~components/cam/src/utils
      
      USE bound_channel_module, only: bound_channel
      USE bound_extra,          only: bound_normal,bound_vert,bound_sync
      USE halo_q,               only: halo_correc_q
      USE halo_vort,            only: vort_comm_pre,vort_comm_post
      USE halo_z,               only: halo_correc_z
      USE vort_3d_module,       only: zeta_diag
      USE wind_module,          only: wind_3d_v2,wind_3d_v1,wind_3d_v0
      USE utils,                only: spline,splint
      USE force_3d_module,      only: add_rperturb

IMPLICIT NONE
PRIVATE

PUBLIC :: init_common,init_channel

CONTAINS

! Local Subroutines:
!----------------------------------------------------------------------------------
! Subroutine init_common           : Initialize common parameters and vertical profiles 
!   Subroutine q3d_vvm_basic_state : Calculate the basic profiles (called by INIT_common) 
!      Subroutine profiles_ussa    : Prepare data (temp. and pressure) 
!      Subroutine vertical_grid    : Determine the vertical grids of VVM
!      Subroutine coords_2d        : Calculate ZZ, ZT, FNZ, FNT
!
! Subroutine init_channel          : Initialize the channel data
!      Subroutine topo_prepare     : Initialize the parameters related to topography  
!----------------------------------------------------------------------------------
 
!===================================================================================
   SUBROUTINE INIT_common
!===================================================================================
      USE mphysics_variables, only: allocate_mphysics_vars
      USE rrtm_vars,          only: allocate_rrtm_vars
      USE rrtm_params,        only: allocate_rrtm_params
      USE trace_gases,        only: allocate_trace_gases 

      ! Local
      INTEGER :: K

      REAL (kind=r8) :: ALPHAW
      REAL (KIND=r8), DIMENSION(NK3) :: TV,ALPHA

      CALL allocate_mphysics_vars()
      CALL allocate_rrtm_vars()
      CALL allocate_rrtm_params()
      CALL allocate_trace_gases()

!***********************************************************************************
!     SET DEFAULT VALUES in CONSTLD (Reset or Recalculated at RESTART)
!***********************************************************************************

!------------------------
! Logical Parameters
!------------------------
      NOSFX      = .FALSE.              ! T, no surface flux
      NOTOPO     = .TRUE.               ! T, no topography
      NOMAP      = .FALSE.              ! T, no (horizontal) mapping
      ZCONST     = .FALSE.              ! T, use of constant-depth vertical grids
      PHYSICS    = .TRUE.               ! T, include physical processes
      TURBULENCE = .FALSE.              ! T, include turbulence
      RADCODE    = .FALSE.              ! T, include radiation (RRTMG)
      MICROCODE  = .TRUE.               ! T, include microphysics
      BUOY       = .TRUE.               ! T, include buoyancy
      NOTHERM    = .FALSE.              ! T, no thermodynamical processes
      FCONST     = .TRUE.               ! T, no Coriolis effect (f=0)

      LOCALP     = .TRUE.               ! T, use local p (& pi) in physics
!------------------------
! REAL Parameters
!------------------------
      DT_GCM  = REAL(get_step_size(), kind=r8)   ! dt of gcm (s) 
      
      DT = 10.0_r8                               ! integration timestep of crm (s)
      
      RXTAU_q = 1800.0_r8      ! relaxation timescale for thermodynamic variables (s) 
      RXTAU_d = 1800.0_r8      ! relaxation timescale for dynamical variables (s)
                               ! If both are specified to zero value, 
                               ! use CAL_RTIME to calculate the relaxation time-scale.
      
!     DZ, DZSQ are determined from SUBROUTINE Q3D_VVM_Basic_State
!     DX, DY, DXSQ, DYSQ, DXDY, DX_GCM, DY_GCM are determined from SUBROUTINE Mapping_channel        
!     WRXMU (=2/DXSQ) is determined from SUBROUTINE Mapping_channel  

      PSFC    = pzero                   ! surface pressure (Pa)
      PISFC   = (PSFC/pzero)**cappa     ! Exner ftn value at sfc (not used so far)
      
      SST     = 298.0_r8                ! sea surface temperature (K)      JUNG: Get this value from CAM-SE
      ZRSEA   = 2.0D-04                 ! surface roughness over sea (m)   JUNG: Get this value from CAM-SE
      ZRLAND  = 1.0D-02                 ! surface roughness over land (m)  JUNG: Get this value from CAM-SE

      ALADV   = 1.0_r8                  ! alpha-value in advection

      DTRAD   = DT_GCM                  ! time interval between radiation calls (s)
      
      CRAD    = 1.0_r8/(3600._r8)       ! Rayleigh-type damping coefficient (s^-1)
!      CRAD    = 0.0_r8            
      
      GWDHT   = 15000.0_r8              ! GWD above this height (m)

      DTHMAX  = 0.5_r8                  ! coefficient for random perturbation (K)
      Z1PERT  = 0.0_r8                  ! z_low for random perturbation  (m)
      Z2PERT  = 100.0_r8                ! z_high for random perturbation (m)

!------------------------
! Integer Parameters
!------------------------
      ITT_MAX = INT(dt_gcm/dt)          ! # of CRM time steps per one GCM time step

      ITINIT  = 0                       ! starting time for random perturbation (currently, not used)
      ITSTOP  = 0                       ! ending time for random perturbation
                                        ! currently: itstop > 0 (add randomperturbation in the initialization) 

!      NSFLUX  = INT(DT_GCM/DT)          ! frequency of surface flux calculation (CRM cycle)
      NSFLUX  = 1                       
      
      NRAD    = INT(DTRAD/DT)           ! frequency of radiation calculation    (CRM cycle)

      NITERW  = 10                      ! # of iterations for w
      NITERXY = 200                     ! # of iterations for psi and chi
      NITERXY_init = 5000               ! # of iterations for psi and chi (initialization)
      
      IVER_wind = 0                     ! Version number of Subroutine wind_3d (0, 1, 2)  

!***********************************************************************************
! Define vertical grids  : ZZ(nk2_ext),ZT(nk2_ext),FNZ(nk2_ext),FNT(nk2_ext)     
! Prepare the basic state: THBAR(nk2_ext)
!                          RHO(nk2_ext),PBAR(nk2_ext),PIBAR(nk2_ext)
!                          RHOZ(nk2_ext),PBARZ(nk2_ext),PIBARZ(nk2_ext)  
!***********************************************************************************

      CALL Q3D_VVM_Basic_State

      DO K=1,NK1
       FN1(K)=RHO(K+1)*FNZ(K)/FNT(K+1)     ! map factor ratio, used in vort_3d_module
       FN2(K)=RHO(K)*FNZ(K)/FNT(K)         ! map factor ratio, used in vort_3d_module
      ENDDO
      FN1(NK2)=RHO(NK2)*FNZ(NK2)/FNT(NK2)  ! not used 
      FN2(NK2)=RHO(NK2)*FNZ(NK2)/FNT(NK2)  ! not used

!     Find the vertical index of the lowest layer where GWD is applied.
      IF (CRAD .NE. 0.0_r8) THEN
        K = 2
        DO WHILE (ZT(K) .LT. GWDHT)
        K = K + 1
        ENDDO
        KLOW_GWD = K
      ENDIF

!***********************************************************************************
!                 Define parameters and vertical profiles used in PHYSICS
!***********************************************************************************

!--------------------------
    IF (physics) THEN
!--------------------------

      ! Initialize Parameters for Microphysics
      CALL initialize_microphysics

      DO K=1,nk1
        RHOL(K) = RHO(K+1)
        PL0(K)  = PBAR(K+1)
        PIL0(K) = PIBAR(K+1)
      ENDDO

      DO K=1,nk1
        DZL(K) = ZZ(K+1) - ZZ(K)
      ENDDO

      RHO_INT(0)   = RHOL(1)
      RHO_INT(nk1) = RHOL(nk1)
      DO K=1,nk1-1
        RHO_INT(K) = &
          (DZL(K)*RHOL(K) + DZL(K+1)*RHOL(K+1))/(2.0_r8*(ZT(K+2) - ZT(K+1)))
      ENDDO

      !--------------------
      IF (radcode) THEN
      !--------------------
      ! Initialize Parameters for Radiation

      ! Get vertical grid parameters and define vertical pressure levels
      PRES(:)  = PBAR(2:nk2_ext) * 0.01_r8    ! [mb]
      PRESI(:) = PBARZ(:) * 0.01_r8           ! [mb]

      NRAD_rad = NRAD

      IF (.NOT.LocalP) &
        CALL TRACE_GAS_INPUT(channel_l, channel_w, nk1_ext, PRES, PRESI, LocalP)
        ! Horizontally uniform profiles are used.

      !--------------------
      ENDIF  ! RADCODE
      !--------------------

!--------------------------
    ENDIF  ! PHYSICS
!--------------------------

   END SUBROUTINE init_common

!===================================================================================
   SUBROUTINE Q3D_VVM_Basic_State
!----------------------------------------------------------------------------------- 
!  INPUT: ZCONST
!  
!  Determine the layer depths: DZ (and DZSQ), DZ_EXT    
!  Determine the vertical grids and map factors: ZZ(nk2_ext),ZT(nk2_ext) 
!                                                FNZ(nk2_ext),FNT(nk2_ext)   
!  Calculate the basic state: THBAR(nk2_ext),
!                             RHO(nk2_ext),PBAR(nk2_ext),PIBAR(nk2_ext)
!                             RHOZ(nk2_ext),PBARZ(nk2_ext),PIBARZ(nk2_ext)  
!
!  nk2_ext = nk1 (28) + num_ext_layers (9) + 1 = 38
!===================================================================================
 
      ! Local parameters and variables
      INTEGER,PARAMETER :: km_hires = 1680
      
      REAL (KIND=r8),DIMENSION(km_hires)  :: height_hires,exner_hires,theta_hires
      REAL (KIND=r8),DIMENSION(nk2_ext) :: thbarz
      
      INTEGER K,L

      ! Compute the high-resolution profiles of exner ftn and theta based on USSAT1976 data. 
      CALL PROFILES_USSA (km_hires,height_hires,exner_hires,theta_hires)

      ! Determine the stretched vertical grids and map factors 
      CALL VERTICAL_GRID (ZCONST,DZ,DZ_EXT,ZZ,ZT,FNZ,FNT)
 
      DZSQ = DZ * DZ
!-----------------------------------------------------------------
!     Determine THBAR, PIBAR, PBAR, RHO at model layer (k = 2~nk2_ext)  
!-----------------------------------------------------------------
      
      do l=2,km_hires
         do k=2,nk2_ext
            if (height_hires(l) .gt. zt(k) .and. height_hires(l-1) .le. zt(k)) then
                THBAR(k) = ((height_hires(l)-zt(k))*theta_hires(l-1)    &
                            +(zt(k)-height_hires(l-1))*theta_hires(l))  &
                           /(height_hires(l)-height_hires(l-1))
            endif
         enddo
      enddo
      
      THBAR(1) = THBAR(2)

      PIBAR(1) = exner_hires(1)
      do l=2,km_hires
         do k=2,nk2_ext
            if (height_hires(l) .gt. zt(k) .and. height_hires(l-1) .le. zt(k)) then
                PIBAR(k) = ((height_hires(l)-zt(k))*exner_hires(l-1)    &
                            +(zt(k)-height_hires(l-1))*exner_hires(l))  &
                           /(height_hires(l)-height_hires(l-1))
            endif
         enddo
      enddo

      do k=1,nk2_ext
         PBAR(k) = pzero*PIBAR(k)**(1.0_r8/cappa)
      enddo

      do k=1,nk2_ext
         RHO(k) = (pzero*PIBAR(k)**((1.0_r8-cappa)/cappa))/(gas_const*THBAR(k))
      enddo

!-----------------------------------------------------------------
!     Determine THBARZ, PIBARZ, PBARZ, RHOZ at model layer (k = 2~nk2)  
!     (THBARZ is not used in the model.  Remove)
!-----------------------------------------------------------------
      
      do l=2,km_hires
         do k=2,nk2_ext
            if (height_hires(l) .gt. zz(k) .and. height_hires(l-1) .le. zz(k)) then
                THBARZ(k) = ((height_hires(l)-zz(k))*theta_hires(l-1)   &
                            +(zz(k)-height_hires(l-1))*theta_hires(l))  &
                           /(height_hires(l)-height_hires(l-1))
            endif
         enddo
      enddo
      
      THBARZ(1) = THBAR(1)
      
      PIBARZ(1) = exner_hires(1)
      do l=2,km_hires
         do k=2,nk2_ext
            if (height_hires(l) .gt. zz(k) .and. height_hires(l-1) .le. zz(k)) then
                PIBARZ(k) = ((height_hires(l)-zz(k))*exner_hires(l-1)   &
                            +(zz(k)-height_hires(l-1))*exner_hires(l))  &
                           /(height_hires(l)-height_hires(l-1))
            endif
         enddo
      enddo

      do k=1,nk2_ext
         PBARZ(k) = pzero*PIBARZ(k)**(1.0_r8/cappa)
      enddo

      do k=1,nk2_ext
         RHOZ(k) = (pzero*PIBARZ(k)**((1.0_r8-cappa)/cappa))/(gas_const*THBARZ(k))
      enddo
            
      CONTAINS

!=======================================================================================      
      SUBROUTINE PROFILES_USSA (km_hires,height_hires,exner_hires,theta_hires)
!=======================================================================================  
!--  This code determines basic state profiles (high resolution) from USSAT1976 data.

      INTEGER,INTENT(IN) :: km_hires
      REAL (KIND=r8),DIMENSION(km_hires),INTENT(OUT) :: height_hires,exner_hires,theta_hires
      
      ! LOCAL
      INTEGER,PARAMETER :: km = 85

      REAL (KIND=r8),DIMENSION(km) :: t_ussatm1976,p_ussatm1976,exner,theta,height, &
                                     temperature,pressure,y2
      
      REAL (KIND=r8),DIMENSION(km_hires) :: temperature_hires,pressure_hires,temp 

      REAL (KIND=r8) :: ratio,yp1,ypkm
      
      INTEGER K

!-- US Standard Atmosphere 1976 pressure and temperature (dz=1000m) upto 85 km ----
      p_ussatm1976(  1 )=   101330.000000000  ;   t_ussatm1976(  1 )=   288.150000000000     
      p_ussatm1976(  2 )=   89880.6900000000  ;   t_ussatm1976(  2 )=   281.651000000000     
      p_ussatm1976(  3 )=   79505.3500000000  ;   t_ussatm1976(  3 )=   275.154100000000    
      p_ussatm1976(  4 )=   70124.6000000000  ;   t_ussatm1976(  4 )=   268.659200000000      
      p_ussatm1976(  5 )=   61663.4600000000  ;   t_ussatm1976(  5 )=   262.166300000000     
      p_ussatm1976(  6 )=   54050.9000000000  ;   t_ussatm1976(  6 )=   255.675500000000     
      p_ussatm1976(  7 )=   47219.9100000000  ;   t_ussatm1976(  7 )=   249.186700000000     
      p_ussatm1976(  8 )=   41107.2200000000  ;   t_ussatm1976(  8 )=   242.700000000000    
      p_ussatm1976(  9 )=   35653.2800000000  ;   t_ussatm1976(  9 )=   236.215200000000     
      p_ussatm1976( 10 )=   30802.1000000000  ;   t_ussatm1976( 10 )=   229.732500000000    
      p_ussatm1976( 11 )=   26501.0800000000  ;   t_ussatm1976( 11 )=   223.251900000000     
      p_ussatm1976( 12 )=   22700.9500000000  ;   t_ussatm1976( 12 )=   216.773300000000     
      p_ussatm1976( 13 )=   19400.2700000000  ;   t_ussatm1976( 13 )=   216.650000000000     
      p_ussatm1976( 14 )=   16580.3100000000  ;   t_ussatm1976( 14 )=   216.650000000000     
      p_ussatm1976( 15 )=   14170.9500000000  ;   t_ussatm1976( 15 )=   216.650000000000     
      p_ussatm1976( 16 )=   12112.2900000000  ;   t_ussatm1976( 16 )=   216.650000000000     
      p_ussatm1976( 17 )=   10353.2200000000  ;   t_ussatm1976( 17 )=   216.650000000000     
      p_ussatm1976( 18 )=   8850.04700000000  ;   t_ussatm1976( 18 )=   216.650000000000     
      p_ussatm1976( 19 )=   7565.49100000000  ;   t_ussatm1976( 19 )=   216.650000000000     
      p_ussatm1976( 20 )=   6467.70300000000  ;   t_ussatm1976( 20 )=   216.650000000000     
      p_ussatm1976( 21 )=   5529.47900000000  ;   t_ussatm1976( 21 )=   216.650000000000     
      p_ussatm1976( 22 )=   4729.07800000000  ;   t_ussatm1976( 22 )=   217.581000000000     
      p_ussatm1976( 23 )=   4047.60600000000  ;   t_ussatm1976( 23 )=   218.574300000000     
      p_ussatm1976( 24 )=   3466.95400000000  ;   t_ussatm1976( 24 )=   219.567200000000     
      p_ussatm1976( 25 )=   2971.81600000000  ;   t_ussatm1976( 25 )=   220.559900000000     
      p_ussatm1976( 26 )=   2549.27500000000  ;   t_ussatm1976( 26 )=   221.552200000000     
      p_ussatm1976( 27 )=   2188.41900000000  ;   t_ussatm1976( 27 )=   222.544300000000     
      p_ussatm1976( 28 )=   1880.00700000000  ;   t_ussatm1976( 28 )=   223.536000000000     
      p_ussatm1976( 29 )=   1616.22000000000  ;   t_ussatm1976( 29 )=   224.527400000000     
      p_ussatm1976( 30 )=   1390.43600000000  ;   t_ussatm1976( 30 )=   225.518600000000     
      p_ussatm1976( 31 )=   1197.04300000000  ;   t_ussatm1976( 31 )=   226.509300000000     
      p_ussatm1976( 32 )=   1031.26800000000  ;   t_ussatm1976( 32 )=   227.499800000000     
      p_ussatm1976( 33 )=   889.067000000000  ;   t_ussatm1976( 33 )=   228.490000000000     
      p_ussatm1976( 34 )=   767.311000000000  ;   t_ussatm1976( 34 )=   230.973700000000     
      p_ussatm1976( 35 )=   663.412000000000  ;   t_ussatm1976( 35 )=   233.744500000000     
      p_ussatm1976( 36 )=   574.593000000000  ;   t_ussatm1976( 36 )=   236.514400000000     
      p_ussatm1976( 37 )=   498.520000000000  ;   t_ussatm1976( 37 )=   239.283400000000     
      p_ussatm1976( 38 )=   433.245000000000  ;   t_ussatm1976( 38 )=   242.051600000000     
      p_ussatm1976( 39 )=   377.135000000000  ;   t_ussatm1976( 39 )=   244.818900000000     
      p_ussatm1976( 40 )=   328.817000000000  ;   t_ussatm1976( 40 )=   247.585400000000    
      p_ussatm1976( 41 )=   287.139000000000  ;   t_ussatm1976( 41 )=   250.351000000000     
      p_ussatm1976( 42 )=   251.128000000000  ;   t_ussatm1976( 42 )=   253.115700000000     
      p_ussatm1976( 43 )=   219.963000000000  ;   t_ussatm1976( 43 )=   255.879600000000     
      p_ussatm1976( 44 )=   192.947000000000  ;   t_ussatm1976( 44 )=   258.642600000000     
      p_ussatm1976( 45 )=   169.492000000000  ;   t_ussatm1976( 45 )=   261.404700000000     
      p_ussatm1976( 46 )=   149.097000000000  ;   t_ussatm1976( 46 )=   264.166000000000     
      p_ussatm1976( 47 )=   131.337000000000  ;   t_ussatm1976( 47 )=   266.926400000000     
      p_ussatm1976( 48 )=   115.847000000000  ;   t_ussatm1976( 48 )=   269.686000000000     
      p_ussatm1976( 49 )=   102.292000000000  ;   t_ussatm1976( 49 )=   270.650000000000     
      p_ussatm1976( 50 )=   90.3330000000000  ;   t_ussatm1976( 50 )=   270.650000000000     
      p_ussatm1976( 51 )=   79.7760000000000  ;   t_ussatm1976( 51 )=   270.650000000000     
      p_ussatm1976( 52 )=   70.4550000000000  ;   t_ussatm1976( 52 )=   270.650000000000     
      p_ussatm1976( 53 )=   62.2120000000000  ;   t_ussatm1976( 53 )=   269.029100000000     
      p_ussatm1976( 54 )=   54.8710000000000  ;   t_ussatm1976( 54 )=   266.274700000000    
      p_ussatm1976( 55 )=   48.3350000000000  ;   t_ussatm1976( 55 )=   263.521200000000     
      p_ussatm1976( 56 )=   42.5220000000000  ;   t_ussatm1976( 56 )=   260.768500000000     
      p_ussatm1976( 57 )=   37.3590000000000  ;   t_ussatm1976( 57 )=   258.016700000000    
      p_ussatm1976( 58 )=   32.7790000000000  ;   t_ussatm1976( 58 )=   255.265700000000     
      p_ussatm1976( 59 )=   28.7210000000000  ;   t_ussatm1976( 59 )=   252.515600000000     
      p_ussatm1976( 60 )=   25.1300000000000  ;   t_ussatm1976( 60 )=   249.766300000000     
      p_ussatm1976( 61 )=   21.9570000000000  ;   t_ussatm1976( 61 )=   247.017900000000     
      p_ussatm1976( 62 )=   19.1560000000000  ;   t_ussatm1976( 62 )=   244.270300000000     
      p_ussatm1976( 63 )=   16.6870000000000  ;   t_ussatm1976( 63 )=   241.523600000000     
      p_ussatm1976( 64 )=   14.5140000000000  ;   t_ussatm1976( 64 )=   238.777800000000     
      p_ussatm1976( 65 )=   12.6040000000000  ;   t_ussatm1976( 65 )=   236.032800000000     
      p_ussatm1976( 66 )=   10.9280000000000  ;   t_ussatm1976( 66 )=   233.288700000000     
      p_ussatm1976( 67 )=   9.46000000000000  ;   t_ussatm1976( 67 )=   230.545400000000     
      p_ussatm1976( 68 )=   8.17400000000000  ;   t_ussatm1976( 68 )=   227.802900000000     
      p_ussatm1976( 69 )=   7.05200000000000  ;   t_ussatm1976( 69 )=   225.061400000000     
      p_ussatm1976( 70 )=   6.07300000000000  ;   t_ussatm1976( 70 )=   222.320600000000     
      p_ussatm1976( 71 )=   5.22000000000000  ;   t_ussatm1976( 71 )=   219.580700000000     
      p_ussatm1976( 72 )=   4.47900000000000  ;   t_ussatm1976( 72 )=   216.841700000000     
      p_ussatm1976( 73 )=   3.83600000000000  ;   t_ussatm1976( 73 )=   214.259700000000     
      p_ussatm1976( 74 )=   3.28000000000000  ;   t_ussatm1976( 74 )=   212.304400000000     
      p_ussatm1976( 75 )=   2.80000000000000  ;   t_ussatm1976( 75 )=   210.349800000000     
      p_ussatm1976( 76 )=   2.38800000000000  ;   t_ussatm1976( 76 )=   208.395800000000     
      p_ussatm1976( 77 )=   2.03300000000000  ;   t_ussatm1976( 77 )=   206.442400000000     
      p_ussatm1976( 78 )=   1.72800000000000  ;   t_ussatm1976( 78 )=   204.489600000000     
      p_ussatm1976( 79 )=   1.46700000000000  ;   t_ussatm1976( 79 )=   202.537400000000     
      p_ussatm1976( 80 )=   1.24300000000000  ;   t_ussatm1976( 80 )=   200.585800000000     
      p_ussatm1976( 81 )=   1.05200000000000  ;   t_ussatm1976( 81 )=   198.634800000000     
      p_ussatm1976( 82 )=  0.889000000000000  ;   t_ussatm1976( 82 )=   196.684400000000     
      p_ussatm1976( 83 )=  0.750000000000000  ;   t_ussatm1976( 83 )=   194.734600000000     
      p_ussatm1976( 84 )=  0.631000000000000  ;   t_ussatm1976( 84 )=   192.785500000000     
      p_ussatm1976( 85 )=  0.531000000000000  ;   t_ussatm1976( 85 )=   190.836900000000 
!----------------------------------------------------------------------------------

      height(1) = 0.0_r8
      do k=2,km
       height(k) = height(k-1) + 1000.0_r8
      enddo

      ratio = psfc/p_ussatm1976(1)       ! psfc is set to pzero

! JUNG_DEBUG (CAM-SE: Pint(sfc) > 1000 mb) extrapolation generates much larger values
! Change the interpolation method from (w.r.t. pi) to (w.r.t. Z)        
!     ratio = 106000./p_ussatm1976(1)  
      
      pressure(1) = p_ussatm1976(1)*ratio
      
      do k=2,km
         pressure(k) = p_ussatm1976(k)*(2.0_r8+ratio)/3.0_r8
      enddo

      do k=1,km
         exner(k) = (pressure(k)/pzero)**cappa
      enddo

      do k=1,km
         theta(k) = 1.01_r8*t_ussatm1976(k)/exner(k)
      enddo

!------- Interpolate USST data to high resolution vertical grid (dz=50m) ------------ 

      height_hires(1) = 0.0_r8
      do k=2,km_hires
         height_hires(k) = height_hires(k-1) + 50.0_r8
      enddo

      yp1  = ( (theta(2)-theta(1))/(height(2)-height(1))             &
               + (theta(3)-theta(2))/(height(3)-height(2)) ) /2.0_r8 

      ypkm = ( (theta(km)-theta(km-1))/(height(km)-height(km-1))     & 
               + (theta(km-1)-theta(km-2))/(height(km-1)-height(km-2)) ) /2.0_r8 

      call spline(height,theta,km,yp1,ypkm,y2)

      theta_hires(1) = theta(1)  
      theta_hires(km_hires) = theta(km)
      do k=2,km_hires-1 
         call splint(height,theta,y2,km,height_hires(k),theta_hires(k))
      enddo

      yp1  = ( (exner(2)-exner(1))/(height(2)-height(1))             &
               + (exner(3)-exner(2))/(height(3)-height(2)) ) /2.0_r8 

      ypkm = ( (exner(km)-exner(km-1))/(height(km)-height(km-1))     & 
               + (exner(km-1)-exner(km-2))/(height(km-1)-height(km-2)) ) /2.0_r8 

      call spline(height,exner,km,yp1,ypkm,y2)

      exner_hires(1) = exner(1)  
      exner_hires(km_hires) = exner(km)
      do k=2,km_hires-1 
         call splint(height,exner,y2,km,height_hires(k),exner_hires(k))
      enddo

      END SUBROUTINE profiles_ussa      
            
!=======================================================================================  
      SUBROUTINE VERTICAL_GRID (Zconstant,DZ,DZ_EXT,zz_total,zt_total,fnz_total,fnt_total)
!======================================================================================= 
      LOGICAL, INTENT(IN) :: Zconstant 
      
      REAL (KIND=r8), INTENT(OUT) :: DZ,DZ_EXT    ! active layer depth & extended layer depth [m] 
      
      REAL (KIND=r8), DIMENSION(nk2_ext), INTENT(OUT) :: ZZ_TOTAL   ! interface heights
      REAL (KIND=r8), DIMENSION(nk2_ext), INTENT(OUT) :: ZT_TOTAL   ! mid-layer heights
      REAL (KIND=r8), DIMENSION(nk2_ext), INTENT(OUT) :: FNZ_TOTAL  ! interface map factor
      REAL (KIND=r8), DIMENSION(nk2_ext), INTENT(OUT) :: FNT_TOTAL  ! mid-layer map factor  

!     Local       
!     Specify the vertical configuration      
      REAL (KIND=r8), PARAMETER :: VIRTUAL_ZTOP = 46000.0_r8  ! ~ 2.3 mb top of the CAM-SE  CAM 6
      REAL (KIND=r8), PARAMETER :: DOMAIN = 30000._r8         ! vertical domain size for active CRM, [m]  
      REAL (KIND=r8), PARAMETER :: ZB  = 0.0_r8               ! surface height, [m]                         
      REAL (KIND=r8), PARAMETER :: DZ1 = 75._r8               ! lowest layer physical depth, [m]            
      
      REAL (KIND=r8) :: CZ1,CZ2                               ! conversion coefficients                                        
      
      INTEGER K
                      
      DZ = DOMAIN/REAL(NK1, kind=r8)
      
      IF (Zconstant) THEN
       CZ2 = 0.0_r8
       CZ1 = 1.0_r8
      ELSE
       CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
       CZ1 = 1.0_r8 - CZ2 * DOMAIN
      ENDIF
      
      CALL COORDS_2D ( NK2, CZ1, CZ2, DZ, ZB, ZZ_total(1:nk2), ZT_total(1:nk2), &
                                             FNZ_total(1:nk2),FNT_total(1:nk2) )

      IF (NUM_EXT_LAYERS .GT. 0) THEN 
       DZ_EXT = (VIRTUAL_ZTOP-ZZ_total(nk2)) / REAL(NUM_EXT_LAYERS, kind=r8)
       
       DO K=nk2+1,nk2_ext
        ZZ_total(K) = ZZ_total(NK2) + DZ_EXT*(K-NK2)
        FNZ_total(K) = 1.0_r8
       ENDDO
      
       DO K=nk2+1,nk2_ext
        ZT_total(K) = ZZ_total(NK2) + DZ_EXT*(REAL(K-NK2, kind=r8) - 0.5_r8) 
        FNT_total(K) = 1.0_r8
       ENDDO 
       
      ELSE
       DZ_EXT = 0.0_r8
      ENDIF
      
      END Subroutine vertical_grid

!======================================================================================= 
      SUBROUTINE COORDS_2D ( KDIM, CZ1, CZ2, DZ, ZB, ZZ, ZT, FNZ, FNT )
!======================================================================================= 
  
      INTEGER,INTENT(IN) :: KDIM
      REAL (KIND=r8), INTENT(IN)  :: CZ1,CZ2,DZ,ZB
      REAL (KIND=r8), DIMENSION(kdim), INTENT(OUT) :: ZZ,ZT,FNZ,FNT
      
      INTEGER :: K

!     Set constant-depth grids
       ZZ(1) = ZB
      DO K = 2, KDIM
       ZZ(K) = ZZ(K-1) + DZ
      ENDDO

      ZT(1) = ZZ(1)
      ZT(2) = ZZ(1) + DZ / 2.0_r8

      DO K = 3, KDIM
       ZT(K) = ZT(K-1) + DZ
      ENDDO
      
!     DEFINE TRANSFORMATION FUNCTIONS ( KLEMP & CHEN ,1982 )
!     AND PHYSICAL COORDINATES

      DO K = 1, KDIM
       FNZ(K) = 1.0_r8 / ( CZ1 + 2.0_r8 * CZ2 * ZZ(K) )
       ZZ(K) = ZZ(K) * ( CZ1 + CZ2 * ZZ(K) )
      ENDDO

      DO K = 1, KDIM
       FNT(K) = 1.0_r8 / ( CZ1 + 2.0_r8 * CZ2 * ZT(K) )
       ZT(K) = ZT(K) * ( CZ1 + CZ2 * ZT(K) )
      ENDDO

      END SUBROUTINE COORDS_2D

END SUBROUTINE q3d_vvm_basic_state   

!===================================================================================
   SUBROUTINE INIT_channel (isRestart,channel)
!===================================================================================
!  Initialize the variables of channel segments
      
      logical, intent(in) :: isRestart  ! .true. iff restart run
      type(channel_t), intent(inout) :: channel  ! channel data

      ! Local
      INTEGER :: I,J,K,nt,num_seg,mi1,mj1,mim,mip,mjm,mjp,mip_c,mjp_c
      
! Calculate various coefficients
!*************************************************************************
      DO num_seg = 1, 4
!************************************************************************* 
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment


      ! Calculate the coefficients used in turbulence and diffusion
      DO J = 1, mj1
      DO I = 1, mi1
      channel%seg(num_seg)%COEFX(I,J) = 1./(channel%seg(num_seg)%RG_T(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY(I,J) = 1./(channel%seg(num_seg)%RG_T(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA(I,J) = 1./(channel%seg(num_seg)%RG_T(I,J)*DXDY*4.0_r8)

      channel%seg(num_seg)%COEFX_K(I,J) = 1./(channel%seg(num_seg)%RG_V(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY_K(I,J) = 1./(channel%seg(num_seg)%RG_V(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA_K(I,J) = 1./(channel%seg(num_seg)%RG_V(I,J)*DXDY*4.0_r8)

      channel%seg(num_seg)%COEFX_E(I,J) = 1./(channel%seg(num_seg)%RG_U(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY_E(I,J) = 1./(channel%seg(num_seg)%RG_U(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA_E(I,J) = 1./(channel%seg(num_seg)%RG_U(I,J)*DXDY*4.0_r8)

      channel%seg(num_seg)%COEFX_Z(I,J) = 1./(channel%seg(num_seg)%RG_Z(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY_Z(I,J) = 1./(channel%seg(num_seg)%RG_Z(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA_Z(I,J) = 1./(channel%seg(num_seg)%RG_Z(I,J)*DXDY*4.0_r8)
      ENDDO
      ENDDO

      ! Calculate the coefficients used in relaxation_2d: 
      ! (Eq. is solved for deviation only, so the gradient terms across the channel are removed.)  
      !-----------------------------
      if (mi1 .gt. mj1) then
      !-----------------------------x-array
      DO J = 1, mj1
       DO I = 1, mi1
          channel%seg(num_seg)%RX2D_C1(I,J) = (channel%seg(num_seg)%RGG_V(1,I,J)         &
                                              +channel%seg(num_seg)%RGG_V(1,I+1,J))/DXSQ 

          channel%seg(num_seg)%RX2D_C2(I,J) = (channel%seg(num_seg)%RGG_U(1,I,J)         &
                                              +channel%seg(num_seg)%RGG_U(1,I-1,J))/DXSQ 
       ENDDO
      ENDDO 
      !-----------------------------  
      else
      !-----------------------------y-array
      DO J = 1, mj1
       DO I = 1, mi1
          channel%seg(num_seg)%RX2D_C1(I,J) = (channel%seg(num_seg)%RGG_U(3,I,J)         &
                                              +channel%seg(num_seg)%RGG_U(3,I,J+1))/DYSQ

          channel%seg(num_seg)%RX2D_C2(I,J) = (channel%seg(num_seg)%RGG_V(3,I,J)         &
                                              +channel%seg(num_seg)%RGG_V(3,I,J-1))/DYSQ
       ENDDO
      ENDDO   
      !-----------------------------         
      endif
      !-----------------------------

      ! Calculate the coefficients used in relaxation_3d
      DO K = 2, NK1
       DO J = 1, mj1
        DO I = 1, mi1
          channel%seg(num_seg)%AGAU_CO(I,J,K-1) = -channel%seg(num_seg)%RG_T(I,J) &
                                                  *FNZ(K)*FNT(K)/(DZSQ*RHO(K))

          channel%seg(num_seg)%BGAU_CO(I,J,K-1) = &
                     (WRXMU+(channel%seg(num_seg)%RGG_U(1,I,J)                    &
                            +channel%seg(num_seg)%RGG_U(1,I-1,J))/DXSQ            &
                           +(channel%seg(num_seg)%RGG_V(3,I,J)                    &
                            +channel%seg(num_seg)%RGG_V(3,I,J-1))/DYSQ)/RHOZ(K)   &
            +channel%seg(num_seg)%RG_T(I,J)*FNZ(K)*(FNT(K+1)/RHO(K+1)+FNT(K)/RHO(K))/DZSQ

          channel%seg(num_seg)%CGAU_CO(I,J,K-1) = -channel%seg(num_seg)%RG_T(I,J) &
                                                  *FNZ(K)*FNT(K+1)/(DZSQ*RHO(K+1))
        ENDDO
       ENDDO
      ENDDO
             
!******************************
      ENDDO   ! num_seg
!******************************

! Prepare characteristics of topography
      CALL TOPO_PREPARE (channel)

!==================================================================================           
      IF (.NOT. isRestart) THEN
!==================================================================================  
      ! Calculate and assign the perturbation to TH3D  
      IF (ITSTOP .GT. 0) CALL ADD_RPERTURB ( DTHMAX,Z1PERT,Z2PERT, channel )
      
! Initialize basic segment data
!*************************************************************************
      DO num_seg = 1, 4
!*************************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment

      mim = channel%seg(num_seg)%mim
      mip = channel%seg(num_seg)%mip

      mjm = channel%seg(num_seg)%mjm
      mjp = channel%seg(num_seg)%mjp

      mip_c = channel%seg(num_seg)%mip_c
      mjp_c = channel%seg(num_seg)%mjp_c

! Surface fluxes: initialization is not really necessary.
      IF (PHYSICS) THEN
      DO J = 1,mj1
       DO I = 1,mi1
        channel%seg(num_seg)%SPREC(I,J) = 0.0_r8

        channel%seg(num_seg)%WTH(I,J) = 0.0_r8
        channel%seg(num_seg)%WQV(I,J) = 0.0_r8
       ENDDO
      ENDDO

      DO J = 1,mjp_c
       DO I = 1,mip_c
        channel%seg(num_seg)%UW(I,J) = 0.0_r8
        channel%seg(num_seg)%WV(I,J) = 0.0_r8

        channel%seg(num_seg)%UW_CON(I,J) = 0.0_r8
        channel%seg(num_seg)%WV_CON(I,J) = 0.0_r8
       ENDDO
      ENDDO
      ENDIF

! Surface information
!-------------------------------------------------------------------
!      GWET is used in the calculation of surface fluxes.
!      TG is used in radiation and surface fluxes.
!      ZROUGH is used in the calculation of turbulence coefficient.
!-------------------------------------------------------------------

      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%TG(I,J) = SST          ! getting from CONSTLD
        channel%seg(num_seg)%ZROUGH(I,J) = ZRSEA    ! getting from CONSTLD
        channel%seg(num_seg)%GWET(I,J) = 1.0_r8
       ENDDO
      ENDDO
!      CALL BOUND_NORMAL  (nhalo,channel,TG=.TRUE.,ZROUGH=.TRUE.,GWET=.TRUE.)
!      CALL BOUND_CHANNEL (nhalo_adv,channel,TG=.TRUE.,ZROUGH=.TRUE.,GWET=.TRUE.)
!      IF (.not.nomap) &
!      CALL HALO_CORREC_Q (nhalo_adv,TG=.TRUE.,ZROUGH=.TRUE.,GWET=.TRUE.)

! Thermodynamic Fields (Temporary setting: zero deviation)

      IF (ITSTOP .GT. 0) THEN
        ! add random perturbation
        DO K = 2, NK2
         DO J = 1, mj1
          DO I = 1, mi1
           channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,K) &
                                            + channel%seg(num_seg)%TH3D_bg(I,J,K)
          ENDDO
         ENDDO
        ENDDO
      ELSE
        DO K = 2, NK2
         DO J = 1, mj1
          DO I = 1, mi1
           channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D_bg(I,J,K)
          ENDDO
         ENDDO
        ENDDO
      ENDIF
      
      DO K = 2, NK2
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%QV3D(I,J,K) = channel%seg(num_seg)%QV3D_bg(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      
      IF (PHYSICS) THEN
      DO K = 2, NK2
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%QC3D(I,J,K) = channel%seg(num_seg)%QC3D_bg(I,J,K)
         channel%seg(num_seg)%QI3D(I,J,K) = channel%seg(num_seg)%QI3D_bg(I,J,K)
         channel%seg(num_seg)%QR3D(I,J,K) = channel%seg(num_seg)%QR3D_bg(I,J,K)
         channel%seg(num_seg)%QS3D(I,J,K) = channel%seg(num_seg)%QS3D_bg(I,J,K)
         channel%seg(num_seg)%QG3D(I,J,K) = channel%seg(num_seg)%QG3D_bg(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      ENDIF

      DO NT = 1, ntracer
      DO K = 2, NK2
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%QT3D(I,J,K,nt) = channel%seg(num_seg)%QT3D_bg(I,J,K,nt)
        ENDDO
       ENDDO
      ENDDO
      ENDDO

!------------------------
      IF (PHYSICS) THEN
!------------------------
      CALL BOUND_NORMAL  (nhalo,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)

      CALL BOUND_CHANNEL (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
      IF (.not.nomap) &
      CALL HALO_CORREC_Q (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)

      CALL BOUND_VERT (channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                       QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
!------------------------
      ELSE
!------------------------
      CALL BOUND_NORMAL  (nhalo,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)
      CALL BOUND_CHANNEL (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)

      IF (.not.nomap) &
      CALL HALO_CORREC_Q (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)

      CALL BOUND_VERT (channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)
!------------------------
      ENDIF  ! PHYSICS
!------------------------

! Vorticity Components (Temporary setting: zero deviation)

      DO K = 2, NK1
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%Z3DX(I,J,K) = channel%seg(num_seg)%Z3DX_bg(I,J,K)
         channel%seg(num_seg)%Z3DY(I,J,K) = channel%seg(num_seg)%Z3DY_bg(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      CALL BOUND_NORMAL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)
      IF (.NOT.nomap) CALL VORT_COMM_PRE (channel)
      CALL BOUND_CHANNEL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)
      IF (.NOT.nomap) CALL VORT_COMM_POST (channel)
      CALL BOUND_VERT (channel,Z3DX=.TRUE.,Z3DY=.TRUE.)

      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ(I,J,nk2) = channel%seg(num_seg)%Z3DZ_bg(I,J,1)
       ENDDO
      ENDDO
      CALL BOUND_NORMAL  (nhalo,channel,Z3DZ=.TRUE.)
      CALL BOUND_SYNC (channel,Z3DZ=.TRUE.)
      CALL BOUND_CHANNEL (nhalo,channel,Z3DZ=.TRUE.)
      IF (.not.nomap) CALL HALO_CORREC_Z (nhalo,channel,Z3DZ=.TRUE.)

      CALL ZETA_DIAG (channel)

! Diagnose wind components

      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%W3D(I,J,  1) = 0.0_r8
        channel%seg(num_seg)%W3D(I,J,nk2) = 0.0_r8
       ENDDO
      ENDDO
      DO K = 2, NK1
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%W3D(I,J,K) = channel%seg(num_seg)%W3D_bg(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      CALL BOUND_NORMAL  (nhalo,channel,W3D=.TRUE.)
      CALL BOUND_CHANNEL (nhalo,channel,W3D=.TRUE.)
      IF (.not.nomap) CALL HALO_CORREC_Q (nhalo,channel,W3D=.TRUE.)

      DO K = 1, NK2
       DO J = mjm, mjp
        DO I = mim, mip
         channel%seg(num_seg)%W3Dnm1(I,J,K) = channel%seg(num_seg)%W3D(I,J,K)
        ENDDO
       ENDDO
      ENDDO

      ! PSI & CHI are only for deviation parts
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%PSI(I,J) = 0.0_r8
         channel%seg(num_seg)%CHI(I,J) = 0.0_r8
       ENDDO
      ENDDO

      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%PSInm1(I,J) = channel%seg(num_seg)%PSI(I,J)
         channel%seg(num_seg)%CHInm1(I,J) = channel%seg(num_seg)%CHI(I,J)
       ENDDO
      ENDDO

!******************************
      ENDDO   ! num_seg
!******************************

#ifdef JUNG_PLAN
      ! UCHANGE=.T., WCHANGE=.T.
!      CALL WIND_3D  (.TRUE., .TRUE., channel)
#else      

      ! UCHANGE=.T., WCHANGE=.T., NITER2D=0: Skip solving the 2D elliptic eq. (deviation = 0)
!      CALL WIND_3D  (.TRUE., .TRUE., channel, NITER2D=0)
      
      SELECT CASE (iver_wind)
      CASE(0)
        CALL WIND_3D_V0  (.TRUE., .TRUE., channel, NITER2D=0)  ! UCHANGE=.T., WCHANGE=.T. 
      CASE(1)
        CALL WIND_3D_V1  (.TRUE., .TRUE., channel, NITER2D=0)  ! UCHANGE=.T., WCHANGE=.T. 
      CASE(2)
        CALL WIND_3D_V2  (.TRUE., .TRUE., channel, NITER2D=0)  ! UCHANGE=.T., WCHANGE=.T. 
      END SELECT      
#endif      

!==================================================================================
      ENDIF   ! .NOT. isRestart 
!==================================================================================  

!*************************************************************************
      DO num_seg = 1, 4
!*************************************************************************
      mip_c = channel%seg(num_seg)%mip_c
      mjp_c = channel%seg(num_seg)%mjp_c

! Surface information
!-------------------------------------------------------------------
!     LOCEAN is used in the calculation of surface fluxes.
!     Formulation can be changed later.
!-------------------------------------------------------------------
      DO J = 1, mjp_c
       DO I = 1, mip_c
        channel%seg(num_seg)%LOCEAN(I,J) = .TRUE.
       ENDDO
      ENDDO
!******************************
      ENDDO   ! num_seg
!******************************

      CONTAINS

!==========================================
      SUBROUTINE TOPO_PREPARE (channel)
!==========================================
!     INPUT : topographic data (height)
!     OUTPUT: associated parameters
!
!     For a while, "NOTOPO=.T." setting is used because there is no input.

      type(channel_t), intent(inout) :: channel  ! channel data

      ! Local variables
      INTEGER :: mi1,mj1,mim,mim_a,mim_c,mip,mip_a,mjm,mjm_a,mjm_c,mjp,mjp_a
      INTEGER :: I,J,K

      INTEGER,DIMENSION(4)    :: KLOWQ_MAX
      INTEGER,DIMENSION(nk2,4):: NTOPOQ,NTOPOU,NTOPOV

      REAL (kind=r8) :: HEIGHT

!*************************************************************************
      IF (NOTOPO) then
!*************************************************************************
!===============================
      DO num_seg = 1, 4
!===============================
      mim   = channel%seg(num_seg)%mim
      mip   = channel%seg(num_seg)%mip
      mim_a = channel%seg(num_seg)%mim_a
      mip_a = channel%seg(num_seg)%mip_a

      mjm   = channel%seg(num_seg)%mjm
      mjp   = channel%seg(num_seg)%mjp
      mjm_a = channel%seg(num_seg)%mjm_a
      mjp_a = channel%seg(num_seg)%mjp_a

      ! HX: level index of topography
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%HX(I,J) = 1
       ENDDO
      ENDDO

      ! KLOWQ: lowest air layer for Q
      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%KLOWQ_IJ(I,J) = channel%seg(num_seg)%HX(I,J)+1
       ENDDO
      ENDDO
      ! KLOWU: lowest air level for U, ETA
      DO J = mjm, mjp
       DO I = mim, mip_a
        channel%seg(num_seg)%KLOWU_IJ(I,J) = channel%seg(num_seg)%HX(I,J)+1
       ENDDO
      ENDDO
      ! KLOWV: lowest air level for V, KSI
      DO J = mjm, mjp_a
       DO I = mim, mip
        channel%seg(num_seg)%KLOWV_IJ(I,J) = channel%seg(num_seg)%HX(I,J)+1
       ENDDO
      ENDDO

      channel%seg(num_seg)%KLOWQ_MAX_GLOB = 2

!-----------------------------------------
! For vorticity advection with mountain
!-----------------------------------------
      channel%seg(num_seg)%DM_KSI_XE = 0
      channel%seg(num_seg)%DM_KSI_XC = 0
      channel%seg(num_seg)%DM_KSI_YS = 0

      channel%seg(num_seg)%DM_ETA_XW = 0
      channel%seg(num_seg)%DM_ETA_YN = 0
      channel%seg(num_seg)%DM_ETA_YC = 0
!--------------------------------
! For turbulence with mountain
!--------------------------------
      channel%seg(num_seg)%DM_Q_TE = 0
      channel%seg(num_seg)%DM_Q_TW = 0
      channel%seg(num_seg)%DM_Q_TN = 0
      channel%seg(num_seg)%DM_Q_TS = 0

      channel%seg(num_seg)%DM_KSI_TE = 0
      channel%seg(num_seg)%DM_KSI_TW = 0
      channel%seg(num_seg)%DM_KSI_TN = 0
      channel%seg(num_seg)%DM_KSI_TS = 0

      channel%seg(num_seg)%DM_ETA_TE = 0
      channel%seg(num_seg)%DM_ETA_TW = 0
      channel%seg(num_seg)%DM_ETA_TN = 0
      channel%seg(num_seg)%DM_ETA_TS = 0

!===============================
      ENDDO   ! num_seg
!===============================

!*************************************************************************
      ELSE   ! NOTOPO
!*************************************************************************
! Prepare information associated with the topographic data

!===============================
      DO num_seg = 1, 4
!===============================
      mim   = channel%seg(num_seg)%mim
      mip   = channel%seg(num_seg)%mip
      mim_a = channel%seg(num_seg)%mim_a
      mip_a = channel%seg(num_seg)%mip_a

      mjm   = channel%seg(num_seg)%mjm
      mjp   = channel%seg(num_seg)%mjp
      mjm_a = channel%seg(num_seg)%mjm_a
      mjp_a = channel%seg(num_seg)%mjp_a

      ! Read Orographic Data (from INPUT): currently, specify 0.0 km
      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%TOPOZ(I,J) = 0.0_r8
       ENDDO
      ENDDO

      ! Halo filling for TOPO
!-----------------------------------------------
!      CALL BOUND_NORMAL  (nhalo,channel,TOPOZ=.TRUE.)
!      CALL BOUND_CHANNEL (nhalo,channel,TOPOZ=.TRUE.)
!
!      IF (.not.nomap) &
!      CALL HALO_CORREC_Q (nhalo,channel,TOPOZ=.TRUE.)
!-----------------------------------------------

      ! HX (i,j): level index of mountain top (at cell center)
      ! HX=1: no mountain   e.g., HX=4: mountain top is at k=4
      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%HX(I,J) = 1
       ENDDO
      ENDDO

      DO J = mjm, mjp
       DO I = mim, mip
        HEIGHT = channel%seg(num_seg)%TOPOZ(I,J)

        IF (HEIGHT.GT.ZZ(1)) THEN
          K = 2
          DO WHILE (HEIGHT.GE.ZZ(K))
           K=K+1
          ENDDO
          channel%seg(num_seg)%HX(I,J) = &
                 (K-1)+NINT((HEIGHT-ZZ(K-1))/(ZZ(K)-ZZ(K-1)))
        ENDIF
       ENDDO
      ENDDO

      ! KLOWQ: lowest air layer for Q
      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%KLOWQ_IJ(I,J) = channel%seg(num_seg)%HX(I,J)+1
       ENDDO
      ENDDO
      ! KLOWU: lowest air level for U, ETA
      DO J = mjm, mjp
       DO I = mim, mip_a
        channel%seg(num_seg)%KLOWU_IJ(I,J) = &
        MAX0((channel%seg(num_seg)%HX(I,J)+1),(channel%seg(num_seg)%HX(I+1,J)+1))
       ENDDO
      ENDDO
      ! KLOWV: lowest air level for V, KSI
      DO J = mjm, mjp_a
       DO I = mim, mip
        channel%seg(num_seg)%KLOWV_IJ(I,J) = &
        MAX0((channel%seg(num_seg)%HX(I,J)+1),(channel%seg(num_seg)%HX(I,J+1)+1))
       ENDDO
      ENDDO

      KLOWQ_MAX(num_seg) = MAXVAL(channel%seg(num_seg)%KLOWQ_IJ(1:mi1,1:mj1))

!===============================
      ENDDO   ! num_seg
!===============================

      channel%seg(:)%KLOWQ_MAX_GLOB = MAXVAL(KLOWQ_MAX(1:4))

!===============================
      DO num_seg = 1, 4
!===============================
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment

      mim_c = channel%seg(num_seg)%mim_c
      mjm_c = channel%seg(num_seg)%mjm_c

!----------------------------------------
! For vorticity advection with mountain
!----------------------------------------
      ! Variables are defined at the vorticity Flux locations
      channel%seg(num_seg)%DM_KSI_XE = 0
      channel%seg(num_seg)%DM_KSI_XC = 0
      channel%seg(num_seg)%DM_KSI_YS = 0

      channel%seg(num_seg)%DM_ETA_XW = 0
      channel%seg(num_seg)%DM_ETA_YN = 0
      channel%seg(num_seg)%DM_ETA_YC = 0

      ! Used in meridional advection for ksi & Zonal advection for eta
      DO J = mjm_c, mj1
       DO I = mim_c, mi1
        DO K = 2,channel%seg(num_seg)%KLOWQ_MAX_GLOB-1
         IF (K.LE.channel%seg(num_seg)%HX(I,J)) channel%seg(num_seg)%DM_KSI_YS(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J)) channel%seg(num_seg)%DM_ETA_XW(I,J,K) = 1
        ENDDO
       ENDDO
      ENDDO

      ! Used in zonal advection for ksi & Meridional advection for eta
      DO J = mim_c, mj1
       DO I = mim_c, mi1
        DO K=2,channel%seg(num_seg)%KLOWQ_MAX_GLOB-1
         IF (K.LE.channel%seg(num_seg)%HX(I,J).OR.K.LE.channel%seg(num_seg)%HX(I,J+1))     &
                                      channel%seg(num_seg)%DM_KSI_XC(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I+2,J).OR.K.LE.channel%seg(num_seg)%HX(I+2,J+1)) &
                                      channel%seg(num_seg)%DM_KSI_XE(I,J,K) = 1

         IF (K.LE.channel%seg(num_seg)%HX(I,J).OR.K.LE.channel%seg(num_seg)%HX(I+1,J))     &
                                      channel%seg(num_seg)%DM_ETA_YC(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J+2).OR.K.LE.channel%seg(num_seg)%HX(I+1,J+2)) &
                                      channel%seg(num_seg)%DM_ETA_YN(I,J,K) = 1
        ENDDO
       ENDDO
      ENDDO

!---------------------------------------------------------------------
! For turbulence with mountain
!---------------------------------------------------------------------
      channel%seg(num_seg)%DM_Q_TE = 0
      channel%seg(num_seg)%DM_Q_TW = 0
      channel%seg(num_seg)%DM_Q_TN = 0
      channel%seg(num_seg)%DM_Q_TS = 0

      channel%seg(num_seg)%DM_KSI_TE = 0
      channel%seg(num_seg)%DM_KSI_TW = 0
      channel%seg(num_seg)%DM_KSI_TN = 0
      channel%seg(num_seg)%DM_KSI_TS = 0

      channel%seg(num_seg)%DM_ETA_TE = 0
      channel%seg(num_seg)%DM_ETA_TW = 0
      channel%seg(num_seg)%DM_ETA_TN = 0
      channel%seg(num_seg)%DM_ETA_TS = 0

      DO J = 1, mj1
       DO I = 1, mi1
        DO K=2,channel%seg(num_seg)%KLOWQ_MAX_GLOB-1
         IF (K.LE.channel%seg(num_seg)%HX(I+1,J)) channel%seg(num_seg)%DM_Q_TE(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I-1,J)) channel%seg(num_seg)%DM_Q_TW(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J+1)) channel%seg(num_seg)%DM_Q_TN(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J-1)) channel%seg(num_seg)%DM_Q_TS(I,J,K) = 1
        ENDDO
       ENDDO
      ENDDO

      DO J = mjm_c, mj1
       DO I = mim_c, mi1
        DO K = 2, channel%seg(num_seg)%KLOWQ_MAX_GLOB-1
         IF (K.LE.channel%seg(num_seg)%HX(I+1,J).OR.K.LE.channel%seg(num_seg)%HX(I+1,J+1)) &
                                      channel%seg(num_seg)%DM_KSI_TE(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I-1,J).OR.K.LE.channel%seg(num_seg)%HX(I-1,J+1)) &
                                      channel%seg(num_seg)%DM_KSI_TW(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J+2)) channel%seg(num_seg)%DM_KSI_TN(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J-1)) channel%seg(num_seg)%DM_KSI_TS(I,J,K) = 1
        ENDDO
       ENDDO
      ENDDO

      DO J = mjm_c, mj1
       DO I = mim_c, mi1
        DO K = 2, channel%seg(num_seg)%KLOWQ_MAX_GLOB-1
         IF (K.LE.channel%seg(num_seg)%HX(I+2,J)) channel%seg(num_seg)%DM_ETA_TE(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I-1,J)) channel%seg(num_seg)%DM_ETA_TW(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J+1).OR.K.LE.channel%seg(num_seg)%HX(I+1,J+1)) &
                                      channel%seg(num_seg)%DM_ETA_TN(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J-1).OR.K.LE.channel%seg(num_seg)%HX(I+1,J-1)) &
                                      channel%seg(num_seg)%DM_ETA_TS(I,J,K) = 1
        ENDDO
       ENDDO
      ENDDO

      DO K = 2, NK2
       NTOPOQ(K,num_seg) = 0
       NTOPOU(K,num_seg) = 0
       NTOPOV(K,num_seg) = 0
       DO J = 1, MJ1
        DO I = 1, MI1
         IF (K.GE.channel%seg(num_seg)%KLOWQ_IJ(I,J)) THEN
          NTOPOQ(K,num_seg) = NTOPOQ(K,num_seg)+1
         ENDIF
         IF (K.GE.channel%seg(num_seg)%KLOWU_IJ(I,J)) THEN
          NTOPOU(K,num_seg) = NTOPOU(K,num_seg)+1
         ENDIF
         IF (K.GE.channel%seg(num_seg)%KLOWV_IJ(I,J)) THEN
          NTOPOV(K,num_seg) = NTOPOV(K,num_seg)+1
         ENDIF
        ENDDO
       ENDDO
      ENDDO

!===============================
      ENDDO   ! num_seg
!===============================

      DO K = 2, NK2
       channel%seg(1)%NTOPOQ_GLOB(K) = NTOPOQ(K,1)+NTOPOQ(K,2)+NTOPOQ(K,3)+NTOPOQ(K,4)
       channel%seg(1)%NTOPOU_GLOB(K) = NTOPOU(K,1)+NTOPOU(K,2)+NTOPOU(K,3)+NTOPOU(K,4)
       channel%seg(1)%NTOPOV_GLOB(K) = NTOPOV(K,1)+NTOPOV(K,2)+NTOPOV(K,3)+NTOPOV(K,4)
      ENDDO

      DO num_seg = 2, 4
       DO K = 2, NK2
        channel%seg(num_seg)%NTOPOQ_GLOB(K) = channel%seg(1)%NTOPOQ_GLOB(K)
        channel%seg(num_seg)%NTOPOU_GLOB(K) = channel%seg(1)%NTOPOU_GLOB(K)
        channel%seg(num_seg)%NTOPOV_GLOB(K) = channel%seg(1)%NTOPOV_GLOB(K)
       ENDDO
      ENDDO

!*************************************************************************
      ENDIF   ! NO_TOPO
!*************************************************************************

      END SUBROUTINE topo_prepare

   END SUBROUTINE init_channel

END MODULE q3d_init_module
