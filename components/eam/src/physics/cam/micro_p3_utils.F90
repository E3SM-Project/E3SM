module micro_p3_utils

  use physics_utils, only: rtype, rtype8, itype, btype

  implicit none
  private
  save

    public :: get_latent_heat, micro_p3_utils_init, &
              avg_diameter, calculate_incloud_mixingratios

    integer, public :: iulog_e3sm
    logical(btype), public :: masterproc_e3sm

    ! Signaling NaN bit pattern that represents a limiter that's turned off.
    integer(itype), parameter :: limiter_off = int(Z'7FF1111111111111', itype)

    real(rtype), public, parameter :: qsmall = 1.e-14_rtype
    real(rtype), public, parameter :: nsmall = 1.e-16_rtype

    real(rtype) :: latent_heat_vapor, latent_heat_sublim, latent_heat_fusion

    real(rtype),public :: rho_1000mb,rho_600mb,ar,br,f1r,f2r,ecr,rho_h2o,kr,kc,aimm,bimm,rin,mi0,nccnst,  &
       eci,eri,bcn,cpw,cons1,cons2,cons3,cons4,cons5,cons6,cons7,         &
       inv_rho_h2o,inv_dropmass,cp,g,rd,rv,ep_2,inv_cp,   &
       thrd,sxth,piov3,piov6,rho_rimeMin,     &
       rho_rimeMax,inv_rho_rimeMax,max_total_ni,dbrk,nmltratio,clbfact_sub,  &
       clbfact_dep, &
       p3_qc_autocon_expon, p3_qc_accret_expon

    real(rtype),dimension(16), public :: dnu

    real(rtype), public, parameter :: mu_r_constant = 1.0_rtype
    real(rtype), public, parameter :: lookup_table_1a_dum1_c =  4.135985029041767d+00 ! 1.0/(0.1*log10(261.7))

    real(rtype),public :: T_zerodegc  ! Temperature at zero degree celcius ~K
    real(rtype),public :: T_rainfrz  ! Contact and immersion freexing temp, -4C  ~K
    real(rtype),public :: T_homogfrz ! Homogeneous freezing temperature, -40C  ~K
    real(rtype),public :: T_icenuc   ! Ice nucleation temperature, -5C ~K

    real(rtype),public :: pi_e3sm
    ! ice microphysics lookup table array dimensions
    integer, public,parameter :: isize        = 50
    integer, public,parameter :: densize      =  5
    integer, public,parameter :: rimsize      =  4
    integer, public,parameter :: rcollsize    = 30
    integer, public,parameter :: ice_table_size      = 12  ! number of quantities used from lookup table
    integer, public,parameter :: collect_table_size  =  2  ! number of ice-rain collection  quantities used from lookup table
    ! switch for warm-rain parameterization
    ! = 1 Seifert and Beheng 2001
    ! = 2 Beheng 1994
    ! = 3 Khairoutdinov and Kogan 2000
    integer, public,parameter :: iparam = 3

    real(rtype), parameter, public :: mincld=0.0001_rtype
    real(rtype), parameter, public :: rho_h2os = 917._rtype  ! bulk density water solid
    real(rtype), parameter, public :: dropmass = 5.2e-7_rtype

    ! particle mass-diameter relationship
    ! currently we assume spherical particles for cloud ice/snow
    ! m = cD^d
    ! exponent
    real(rtype), parameter :: dsph = 3._rtype

    ! Bounds for mean diameter for different constituents.
    ! real(rtype), parameter :: lam_bnd_rain(2) = 1._rtype/[500.e-6_rtype, 20.e-6_rtype]
    ! real(rtype), parameter :: lam_bnd_snow(2) = 1._rtype/[2000.e-6_rtype, 10.e-6_rtype]

    ! Minimum average mass of particles.
    real(rtype), parameter :: min_mean_mass_liq = 1.e-20_rtype
    real(rtype), parameter :: min_mean_mass_ice = 1.e-20_rtype

    ! in-cloud values
    REAL(rtype), PARAMETER :: min_cld_frac   = 1.e-20_rtype !! threshold min value for cloud fraction
    real(rtype), parameter :: incloud_limit = 5.1E-3
    real(rtype), parameter :: precip_limit  = 1.0E-2

    contains
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine micro_p3_utils_init(cpair,rair,rh2o,rhoh2o,mwh2o,mwdry,gravit,latvap,latice, &
                   cpliq,tmelt,pi,iulog,masterproc)

    real(rtype), intent(in) :: cpair
    real(rtype), intent(in) :: rair
    real(rtype), intent(in) :: rh2o
    real(rtype), intent(in) :: rhoh2o
    real(rtype), intent(in) :: mwh2o
    real(rtype), intent(in) :: mwdry
    real(rtype), intent(in) :: gravit
    real(rtype), intent(in) :: latvap
    real(rtype), intent(in) :: latice
    real(rtype), intent(in) :: cpliq
    real(rtype), intent(in) :: tmelt
    real(rtype), intent(in) :: pi
    integer, intent(in)     :: iulog
    logical(btype), intent(in)     :: masterproc

    ! logfile info
    iulog_e3sm      = iulog
    masterproc_e3sm = masterproc

    ! mathematical/optimization constants
    thrd  = 1._rtype/3._rtype
    sxth  = 1._rtype/6._rtype 
    pi_e3sm = pi
    piov3 = pi*thrd
    piov6 = pi*sxth

    ! maximum total ice concentration (sum of all categories)
     max_total_ni = 500.e+3_rtype  !(m)

    ! droplet concentration (m-3)
    nccnst = 200.e+6_rtype

    ! parameters for Seifert and Beheng (2001) autoconversion/accretion
    kc     = 9.44e+9_rtype
    kr     = 5.78e+3_rtype

    ! Temperature parameters
    T_zerodegc  = tmelt 
    T_homogfrz = tmelt-40._rtype
    T_icenuc   = tmelt-15._rtype
    T_rainfrz  = tmelt-4._rtype

    ! physical constants
    cp     = cpair ! specific heat of dry air (J/K/kg) !1005.
    inv_cp = 1._rtype/cp ! inverse of cp
    g      = gravit ! Gravity (m/s^2) !9.816
    rd     = rair ! Dry air gas constant     ~ J/K/kg
    rv     = rh2o ! Water vapor gas constant ~ J/K/kg     !461.51
    ep_2   = mwh2o/mwdry  ! ratio of molecular mass of water to the molecular mass of dry air !0.622
    rho_1000mb = 100000._rtype/(rd*T_zerodegc) ! density of air at surface
    rho_600mb = 60000._rtype/(rd*253.15_rtype)
    ar     = 841.99667_rtype 
    br     = 0.8_rtype
    f1r    = 0.78_rtype
    f2r    = 0.32_rtype
    ecr    = 1._rtype
    rho_h2o   = rhoh2o ! Density of liquid water (STP) !997.
    cpw    = cpliq  ! specific heat of fresh h2o (J/K/kg) !4218.
    inv_rho_h2o = 1._rtype/rho_h2o  !inverse of (max.) density of liquid water
    inv_dropmass = 1._rtype/dropmass  !inverse of dropmass

    latent_heat_vapor = latvap           ! latent heat of vaporization
    latent_heat_sublim = latvap + latice  ! latent heat of sublimation
    latent_heat_fusion  = latice           ! latent heat of fusion

    ! limits for rime density [kg m-3]
    rho_rimeMin     =  50._rtype
    rho_rimeMax     = 900._rtype
    inv_rho_rimeMax =   1._rtype/rho_rimeMax

    ! Bigg (1953)
    !bimm   = 100.
    !aimm   = 0.66
    ! Barklie and Gokhale (1959)
    bimm   = 2._rtype
    aimm   = 0.65_rtype
    rin    = 0.1e-6_rtype
    mi0    = 4._rtype*piov3*900._rtype*1.e-18_rtype

    eci    = 0.5_rtype
    eri    = 1._rtype
    bcn    = 2._rtype

    ! mean size for soft lambda_r limiter [microns]
    dbrk   = 600.e-6_rtype
    ! ratio of rain number produced to ice number loss from melting
    nmltratio = 1.0_rtype

    cons1 = piov6*rho_h2o
    cons2 = 4._rtype*piov3*rho_h2o
    cons3 = 1._rtype/(cons2*1.562500000000000d-14)  ! 1._rtype/(cons2*bfb_pow(25.e-6_rtype,3.0_rtype))
    cons4 = 1._rtype/(dbrk**3*pi*rho_h2o)
    cons5 = piov6*bimm
    cons6 = piov6**2*rho_h2o*bimm
    cons7 = 4._rtype*piov3*rho_h2o*1.e-18_rtype

    ! droplet spectral shape parameter for mass spectra, used for Seifert and Beheng (2001)
    ! warm rain autoconversion/accretion option only (iparam = 1)
!    allocate(dnu(16))
    dnu(1)  =  0.000_rtype
    dnu(2)  = -0.557_rtype
    dnu(3)  = -0.430_rtype
    dnu(4)  = -0.307_rtype
    dnu(5)  = -0.186_rtype
    dnu(6)  = -0.067_rtype
    dnu(7)  = -0.050_rtype
    dnu(8)  = -0.167_rtype
    dnu(9)  = -0.282_rtype
    dnu(10) = -0.397_rtype
    dnu(11) = -0.512_rtype
    dnu(12) = -0.626_rtype
    dnu(13) = -0.739_rtype
    dnu(14) = -0.853_rtype
    dnu(15) = -0.966_rtype
    dnu(16) = -0.966_rtype

    ! calibration factors for ice deposition and sublimation
    !   These are adjustable ad hoc factors used to increase or decrease deposition and/or
    !   sublimation rates.  The representation of the ice capacitances are highly simplified
    !   and the appropriate values in the diffusional growth equation are uncertain.
    clbfact_dep = 1._rtype
    clbfact_sub = 1._rtype

    return
    end subroutine micro_p3_utils_init
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine get_latent_heat(its,ite,kts,kte,v,s,f)

       integer,intent(in) :: its,ite,kts,kte
       real(rtype),dimension(its:ite,kts:kte),intent(out) :: v,s,f

!       integer i,k

       v(:,:) = latent_heat_vapor !latvap           ! latent heat of vaporization
       s(:,:) = latent_heat_sublim !latvap + latice  ! latent heat of sublimation
       f(:,:) = latent_heat_fusion  !latice           ! latent heat of fusion
 
! Original P3 definition of latent heats:   
!       do i = its,ite
!          do k = kts,kte
!          latent_heat_vapor(i,k)    = 3.1484e6-2370.*t(i,k)
!          latent_heat_sublim(i,k)    = latent_heat_vapor(i,k)+0.3337e6
!          latent_heat_fusion(i,k)     = latent_heat_sublim(i,k)-latent_heat_vapor(i,k)
!          end do
!       end do
       return
    end subroutine get_latent_heat

!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine calculate_incloud_mixingratios(qc,qr,qi,qm,nc,nr,ni,bm, &
          inv_cld_frac_l,inv_cld_frac_i,inv_cld_frac_r, &
          qc_incld,qr_incld,qi_incld,qm_incld,nc_incld,nr_incld,ni_incld,bm_incld)

       real(rtype),intent(in)   :: qc, qr, qi, qm
       real(rtype),intent(in)   :: nc, nr, ni, bm
       real(rtype),intent(in)   :: inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r
       real(rtype),intent(out)  :: qc_incld, qr_incld, qi_incld, qm_incld
       real(rtype),intent(out)  :: nc_incld, nr_incld, ni_incld, bm_incld

       if (qc.ge.qsmall) then
          qc_incld = qc*inv_cld_frac_l
          nc_incld = max(nc*inv_cld_frac_l,0._rtype)
          !AaronDonahue, kai has something about if nccons then nc=ncnst/rho
       else
          qc_incld = 0._rtype
          nc_incld = 0._rtype
       end if 
       if (qi.ge.qsmall) then
          qi_incld = qi*inv_cld_frac_i
          ni_incld = max(ni*inv_cld_frac_i,0._rtype)
          !AaronDonahue, kai has something about if nicons then ni=ninst/rho
       else
          qi_incld = 0._rtype
          ni_incld = 0._rtype
       end if 
       if (qm.ge.qsmall.and.qi.ge.qsmall) then
          qm_incld = qm*inv_cld_frac_i
          bm_incld = max(bm*inv_cld_frac_l,0._rtype)
       else
          qm_incld = 0._rtype
          bm_incld = 0._rtype
       end if 
       if (qr.ge.qsmall) then
          qr_incld = qr*inv_cld_frac_r
          nr_incld = max(nr*inv_cld_frac_r,0._rtype)
          !AaronDonahue, kai has something about if nccons then nc=ncnst/rho
       else
          qr_incld = 0._rtype
          nr_incld = 0._rtype
       end if
       if (qc_incld.gt.incloud_limit .or.qi_incld.gt.incloud_limit &
            .or. qr_incld.gt.precip_limit .or.bm_incld.gt.incloud_limit) then
          !write(errmsg,'(a3,i4,3(a5,1x,e16.8,1x))') 'k: ', k, ', qc:',qc_incld, &
          !     ', qi:',qi_incld,', qr:',qr_incld
          qc_incld    = min(qc_incld,incloud_limit)
          qi_incld = min(qi_incld,incloud_limit)
          bm_incld = min(bm_incld,incloud_limit)
          qr_incld    = min(qr_incld,precip_limit)
!          if (masterproc) write(iulog,*)  errmsg

!          call handle_errmsg('Micro-P3 (Init)',subname='In-cloud mixing
!          ratio too large',extra_msg=errmsg)
       end if
    end subroutine calculate_incloud_mixingratios
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!

real(rtype) elemental function avg_diameter(q, n, rho_air, rho_sub)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  ! Assumes that diameter follows an exponential distribution.
  real(rtype), intent(in) :: q         ! mass mixing ratio
  real(rtype), intent(in) :: n         ! number concentration (per volume)
  real(rtype), intent(in) :: rho_air   ! local density of the air
  real(rtype), intent(in) :: rho_sub   ! density of the particle substance

  avg_diameter = (pi_e3sm * rho_sub * n/(q*rho_air))**(-1._rtype/3._rtype)

end function avg_diameter


end module micro_p3_utils
