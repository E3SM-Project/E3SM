  module radar_simulator_types

! Collection of common variables and types
! Part of QuickBeam v1.03 by John Haynes
! Updated by Roj Marchand June 2010 

  integer, parameter ::       &
  maxhclass = 20         ,& ! max number of hydrometeor classes
  nRe_types = 550       ! max number or Re size bins allowed in N and Z_scaled look up table

  ! These variables define discrete diameters used to represent the DSDs.
  integer, parameter ::       &
  nd = 85               ! number of discrete particles used in construction DSDs
  real*8, parameter ::        &
  dmin = 0.1                 ,& ! min size of discrete particle
  dmax = 10000.                 ! max size of discrete particle
   
  integer, parameter :: &   ! These parameters used to define temperature intervals in mie LUTs
  mt_nfreq = 5              , &
  mt_ntt = 39               , & ! num temperatures in table
  mt_nf = 14            , & ! number of ice fractions in table  
  mt_nd = 85                    ! num discrete mode-p drop sizes in table

  integer, parameter :: &   ! These parameters used to defines Re intervals in scale LUTs
  Re_BIN_LENGTH=10, &        
  Re_MAX_BIN=250

  integer, parameter :: &   ! These parameters used to define Temperature invervals in scale LUTs
  cnt_liq = 19, &       ! liquid temperature count
  cnt_ice = 20          ! ice temperature count


! ---- hydrometeor class type -----  
  
  type class_param
  
    ! variables used to store hydrometeor "default" properties
    real*8,  dimension(maxhclass) :: p1,p2,p3,dmin,dmax,apm,bpm,rho
    integer, dimension(maxhclass) :: dtype,col,cp,phase
  
    ! Radar properties
    real*8  :: freq,k2
    integer :: nhclass      ! number of hydrometeor classes in use
    integer :: use_gas_abs, do_ray
    
    ! defines location of radar relative to hgt_matrix.   
    logical :: radar_at_layer_one	! if true radar is assume to be at the edge 
    					! of the first layer, if the first layer is the
    					! surface than a ground-based radar.   If the
    					! first layer is the top-of-atmosphere, then
    					! a space borne radar. 
    
    ! variables used to store Z scale factors
    character*240 :: scale_LUT_file_name
    logical :: load_scale_LUTs, update_scale_LUTs
    logical, dimension(maxhclass,nRe_types) :: N_scale_flag
    logical, dimension(maxhclass,mt_ntt,nRe_types) :: Z_scale_flag,Z_scale_added_flag
    real*8,  dimension(maxhclass,mt_ntt,nRe_types) :: Ze_scaled,Zr_scaled,kr_scaled
    real*8,  dimension(maxhclass,nd,nRe_types) :: fc, rho_eff

    ! used to determine Re index
    real*8  :: step_list(Re_MAX_BIN),base_list(Re_MAX_BIN)
  
    ! used to determine temperature index
    real*8 :: &
        mt_ttl(cnt_liq), &  ! liquid temperatures (K)
        mt_tti(cnt_ice)     ! ice temperatures (K)

    real*8 :: D(nd) ! set of discrete diameters used to represent DSDs

  end type class_param
  
   
  end module radar_simulator_types
