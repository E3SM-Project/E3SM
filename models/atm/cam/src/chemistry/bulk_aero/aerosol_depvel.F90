!--------------------------------------------------------------------------------
! used for dust and seasalt 
!--------------------------------------------------------------------------------
module aerosol_depvel

contains

!--------------------------------------------------------------------------------
! settling velocity
!--------------------------------------------------------------------------------
subroutine aerosol_depvel_compute( ncol, nlev, naero, t, pmid, ram1, fv, diam, stk_crc, dns_aer, &
                                   vlc_dry, vlc_trb, vlc_grv )

  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst,    only: pi, gravit, rair, boltz

  ! !ARGUMENTS:
  !
  implicit none
  !
  integer,  intent(in) :: ncol,nlev
  integer,  intent(in) :: naero
  real(r8), intent(in) :: t(:,:)          !atm temperature (K)
  real(r8), intent(in) :: pmid(:,:)       !atm pressure (Pa)
  real(r8), intent(in) :: fv(:)           !friction velocity (m/s)
  real(r8), intent(in) :: ram1(:)         !aerodynamical resistance (s/m)
  real(r8), intent(in) :: diam(:,:,:)
  real(r8), intent(in) :: stk_crc(:)
  real(r8), intent(in) :: dns_aer

  real(r8), intent(out) :: vlc_trb(:,:)    !Turbulent deposn velocity (m/s)
  real(r8), intent(out) :: vlc_grv(:,:,:)  !grav deposn velocity (m/s)
  real(r8), intent(out) :: vlc_dry(:,:,:)  !dry deposn velocity (m/s)

  !------------------------------------------------------------------------
  ! Local Variables
  integer  :: m,i,k          !indices
  real(r8) :: vsc_dyn_atm(ncol,nlev)   ![kg m-1 s-1] Dynamic viscosity of air
  real(r8) :: vsc_knm_atm(ncol,nlev)   ![m2 s-1] Kinematic viscosity of atmosphere
  real(r8) :: shm_nbr_xpn   ![frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
  real(r8) :: shm_nbr       ![frc] Schmidt number
  real(r8) :: stk_nbr       ![frc] Stokes number
  real(r8) :: mfp_atm(ncol,nlev)       ![m] Mean free path of air
  real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
  real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
  real(r8) :: slp_crc(ncol,nlev,naero) ![frc] Slip correction factor
  real(r8) :: rss_lmn(naero) ![s m-1] Quasi-laminar layer resistance
  real(r8) :: tmp          !temporary 

  ! constants
  real(r8),parameter::shm_nbr_xpn_lnd=-2._r8/3._r8 ![frc] shm_nbr_xpn over land
  real(r8),parameter::shm_nbr_xpn_ocn=-1._r8/2._r8 ![frc] shm_nbr_xpn over ccean

  real(r8) :: rho                 !atm density (kg/m**3)

  ! needs fv and ram1 passed in from lnd model

  !------------------------------------------------------------------------

  do k=1,nlev
     do i=1,ncol
        rho = pmid(i,k)/rair/t(i,k)
        ! from subroutine dst_dps_dry (consider adding sanity checks from line 212)
        ! when code asks to use midlayer density, pressure, temperature,
        ! I use the data coming in from the atmosphere, ie t(i,k), pmid(i,k)

        ! Quasi-laminar layer resistance: call rss_lmn_get
        ! Size-independent thermokinetic properties
        vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
             (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
        mfp_atm(i,k) = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
             (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))
        vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho ![m2 s-1] Kinematic viscosity of air

        do m = 1, naero
           slp_crc(i,k,m) = 1.0_r8 + 2.0_r8 * mfp_atm(i,k) * &
                (1.257_r8+0.4_r8*exp(-1.1_r8*diam(i,k,m)/(2.0_r8*mfp_atm(i,k)))) / &
                diam(i,k,m)   ![frc] Slip correction factor SeP97 p. 464
           vlc_grv(i,k,m) = (1.0_r8/18.0_r8) * diam(i,k,m) * diam(i,k,m) * dns_aer * &
                gravit * slp_crc(i,k,m) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466
           vlc_grv(i,k,m) = vlc_grv(i,k,m) * stk_crc(m)         ![m s-1] Correction to Stokes settling velocity
           vlc_dry(i,k,m)=vlc_grv(i,k,m)
        end do

     enddo
  enddo
  k=nlev  ! only look at bottom level for next part
  do m = 1, naero
     do i=1,ncol
        stk_nbr = vlc_grv(i,k,m) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))    ![frc] SeP97 p.965
        dff_aer = boltz * t(i,k) * slp_crc(i,k,m) / &    ![m2 s-1]
             (3.0_r8*pi*vsc_dyn_atm(i,k)*diam(i,k,m)) !SeP97 p.474
        shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972
        shm_nbr_xpn = shm_nbr_xpn_lnd                          ![frc]
        !           if(ocnfrac.gt.0.5) shm_nbr_xpn=shm_nbr_xpn_ocn
        ! fxm: Turning this on dramatically reduces
        ! deposition velocity in low wind regimes
        ! Schmidt number exponent is -2/3 over solid surfaces and
        ! -1/2 over liquid surfaces SlS80 p. 1014
        ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd
        ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt # 
        tmp = shm_nbr**shm_nbr_xpn + 10.0_r8**(-3.0_r8/stk_nbr)
        rss_lmn(m) = 1.0_r8 / (tmp*fv(i)) ![s m-1] SeP97 p.972,965

        rss_trb = ram1(i) + rss_lmn(m) + ram1(i)*rss_lmn(m)*vlc_grv(i,k,m) ![s m-1]
        vlc_trb(i,m) = 1.0_r8 / rss_trb                            ![m s-1]
        vlc_dry(i,k,m) = vlc_trb(i,m)  +vlc_grv(i,k,m)
     end do !ncol
  end do

end subroutine aerosol_depvel_compute

end module aerosol_depvel
