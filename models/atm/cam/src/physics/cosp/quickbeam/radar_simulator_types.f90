  module radar_simulator_types

! Collection of common variables and types
! Part of QuickBeam v1.03 by John Haynes
! http://reef.atmos.colostate.edu/haynes/radarsim

  integer, parameter ::       &
  maxhclass = 20         ,& ! max number of hydrometeor classes
  nd = 85           ,& ! number of discrete particles  
  nRe_types = 250        ! number or Re size bins allowed in N and Z_scaled look up table

  real*8, parameter ::        &
  dmin = 0.1                 ,& ! min size of discrete particle
  dmax = 10000.                    ! max size of discrete particle
   
  integer, parameter :: &
  mt_nfreq = 5              , &
  mt_ntt = 39               , &    ! num temperatures in table
  mt_nf   = 14          , &   ! number of ice fractions in table  
  mt_nd = 85                   ! num discrete mode-p drop sizes in table


! ---- hydrometeor class type -----  
  
  type class_param
    real*8,  dimension(maxhclass) :: p1,p2,p3,dmin,dmax,apm,bpm,rho
    integer, dimension(maxhclass) :: dtype,col,cp,phase
    logical, dimension(maxhclass,nRe_types) :: scaled
    logical, dimension(maxhclass,mt_ntt,nRe_types) :: z_flag
    real*8,  dimension(maxhclass,mt_ntt,nRe_types) :: Ze_scaled,Zr_scaled,kr_scaled
    real*8,  dimension(maxhclass,nd,nRe_types) :: fc, rho_eff
    integer, dimension(maxhclass,nd,nRe_types) :: ifc
    integer, dimension(maxhclass) :: idd
  end type class_param

! ----- mie table structure -----
  
  type mie
    real*8 :: freq(mt_nfreq), tt(mt_ntt), f(mt_nf), D(mt_nd)
    real*8, dimension(mt_nd,mt_ntt,mt_nf,mt_nfreq) :: qext, qbsca
    integer :: phase(mt_ntt)
  end type mie

  real*8, dimension(:), allocatable :: &
    mt_qext, mt_qbsca         ! extincion/backscatter efficiency
  
  real*8 :: &
    mt_ttl(20), &        ! liquid temperatures (C)
    mt_tti(19)           ! ice temperatures (C)

  integer*4 :: &
    cnt_liq, &           ! liquid temperature count
    cnt_ice              ! ice temperature count

  end module radar_simulator_types
