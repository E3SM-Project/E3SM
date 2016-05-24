module DUSTMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: DUSTMod
!
! !DESCRIPTION: 
! Routines in this module calculate Dust mobilization and dry deposition for dust.
! Simulates dust mobilization due to wind from the surface into the 
! lowest atmospheric layer. On output flx_mss_vrt_dst(ndst) is the surface dust 
! emission (kg/m**2/s) [ + = to atm].
! Calculates the turbulent component of dust dry deposition, (the turbulent deposition 
! velocity through the lowest atmospheric layer). CAM will calculate the settling 
! velocity through the whole atmospheric column. The two calculations will determine 
! the dust dry deposition flux to the surface.
!                              
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8 
  use clmtype
  use clm_varpar  , only : dst_src_nbr, ndst, sz_nbr
  use clm_varcon  , only : grav, istsoil
  use clm_varcon  , only : istcrop, istice_mec
  use clm_varctl  , only : iulog
  use abortutils  , only : endrun
  use subgridAveMod, only: p2l_1d
  use clm_varcon, only: spval
!  
! !PUBLIC TYPES
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public Dustini        ! Initialize variables used in subroutine Dust
  public DustEmission   ! Dust mobilization 
  public DustDryDep     ! Turbulent dry deposition for dust
!
! !REVISION HISTORY
! Created by Sam Levis, updated to clm2.1 by Mariana Vertenstein
! Source: C. Zender's dust model
!
!EOP
!
! Data private to this module
!
  private
  real(r8) ovr_src_snk_mss(dst_src_nbr,ndst)  
  real(r8) tmp1          !Factor in saltation computation (named as in Charlie's code)
  real(r8) dmt_vwr(ndst) ![m] Mass-weighted mean diameter resolved
  real(r8) stk_crc(ndst) ![frc] Correction to Stokes settling velocity
  real(r8) dns_aer       ![kg m-3] Aerosol density
!------------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DustEmission
!
! !INTERFACE:
  subroutine DustEmission (lbp, ubp, lbc,ubc,lbl,ubl,num_nolakep, filter_nolakep)
!
! !DESCRIPTION: 
! Dust mobilization. This code simulates dust mobilization due to wind
! from the surface into the lowest atmospheric layer
! On output flx_mss_vrt_dst(ndst) is the surface dust emission 
! (kg/m**2/s) [ + = to atm]
! Source: C. Zender's dust model
!
! !USES
   use clm_atmlnd   , only : clm_a2l
   use shr_const_mod, only : SHR_CONST_RHOFW
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp,lbc,ubc,ubl,lbl                    ! pft bounds
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(num_nolakep) ! pft filter for non-lake points
!
! !LOCAL VARIABLES
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)         ! pft's column index
    integer , pointer :: plandunit(:)       ! pft's landunit index
    integer , pointer :: pgridcell(:)       ! pft's gridcell index
    integer , pointer :: ityplun(:)         ! landunit type
    real(r8), pointer :: tlai(:)            ! one-sided leaf area index, no burying by snow
    real(r8), pointer :: tsai(:)            ! one-sided stem area index, no burying by snow
    real(r8), pointer :: frac_sno(:)        ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: gwc_thr(:)         ! threshold gravimetric soil moisture based on clay content
    real(r8), pointer :: forc_rho(:)        ! density (kg/m**3)
    real(r8), pointer :: fv(:)              ! friction velocity (m/s) (for dust model)
    real(r8), pointer :: u10(:)             ! 10-m wind (m/s) (created for dust model)
    real(r8), pointer :: mbl_bsn_fct(:)     ! basin factor
    real(r8), pointer :: mss_frc_cly_vld(:) ! [frc] Mass fraction clay limited to 0.20
    real(r8), pointer :: h2osoi_vol(:,:)    ! volumetric soil water (0<=h2osoi_vol<=watsat)
    real(r8), pointer :: h2osoi_liq(:,:)    ! liquid soil water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)    ! frozen soil water (kg/m2)
    real(r8), pointer :: watsat(:,:)        ! saturated volumetric soil water

! local pointers to implicit out arguments
!
    real(r8), pointer :: flx_mss_vrt_dst(:,:)   ! surface dust emission (kg/m**2/s) 
    real(r8), pointer :: flx_mss_vrt_dst_tot(:) ! total dust flux into atmosphere 

! !REVISION HISTORY
! Created by Sam Levis
! Migrated to new data structures by Peter Thornton and Mariana Vertenstein
! !Created by Peter Thornton and Mariana Vertenstein
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: fp,p,c,l,g,m,n      ! indices
    real(r8) :: liqfrac             ! fraction of total water that is liquid
    real(r8) :: wnd_frc_rat         ! [frc] Wind friction threshold over wind friction
    real(r8) :: wnd_frc_slt_dlt     ! [m s-1] Friction velocity increase from saltatn
    real(r8) :: wnd_rfr_dlt         ! [m s-1] Reference windspeed excess over threshld
    real(r8) :: dst_slt_flx_rat_ttl
    real(r8) :: flx_mss_hrz_slt_ttl
    real(r8) :: flx_mss_vrt_dst_ttl(lbp:ubp)
    real(r8) :: frc_thr_wet_fct
    real(r8) :: frc_thr_rgh_fct
    real(r8) :: wnd_frc_thr_slt
    real(r8) :: wnd_rfr_thr_slt
    real(r8) :: wnd_frc_slt
    real(r8) :: lnd_frc_mbl(lbp:ubp)
    real(r8) :: bd
    real(r8) :: gwc_sfc
    real(r8) :: ttlai(lbp:ubp)
    real(r8) :: tlai_lu(lbl:ubl)
!    
! constants
!
    real(r8), parameter :: cst_slt = 2.61_r8           ! [frc] Saltation constant
    real(r8), parameter :: flx_mss_fdg_fct = 5.0e-4_r8 ! [frc] Empir. mass flx tuning eflx_lh_vegt
    real(r8), parameter :: vai_mbl_thr = 0.3_r8        ! [m2 m-2] VAI threshold quenching dust mobilization
    real(r8), pointer :: wtlunit(:)         ! weight of pft relative to landunit
    real(r8) :: sumwt(lbl:ubl)              ! sum of weights
    logical  :: found                       ! temporary for error check
    integer  :: index

!------------------------------------------------------------------------

    ! Assign local pointers to derived type scalar members (gridcell-level)

    forc_rho        => clm_a2l%forc_rho

    ! Assign local pointers to derived type scalar members (landunit-level)

    ityplun         => lun%itype

    ! Assign local pointers to derived type scalar members (column-level)

    frac_sno        => cps%frac_sno
    gwc_thr         => cps%gwc_thr
    mbl_bsn_fct     => cps%mbl_bsn_fct
    mss_frc_cly_vld => cps%mss_frc_cly_vld
    h2osoi_vol      => cws%h2osoi_vol
    h2osoi_liq      => cws%h2osoi_liq
    h2osoi_ice      => cws%h2osoi_ice
    watsat          => cps%watsat

    ! Assign local pointers to derived type scalar members (pft-level)

    pgridcell       => pft%gridcell
    plandunit       => pft%landunit
    pcolumn         => pft%column
    tlai            => pps%tlai
    tsai            => pps%tsai
    fv              => pps%fv
    u10             => pps%u10
    flx_mss_vrt_dst => pdf%flx_mss_vrt_dst
    flx_mss_vrt_dst_tot => pdf%flx_mss_vrt_dst_tot
   !local pointers from subgridAveMod/p2l_1d
    wtlunit         => pft%wtlunit

    ttlai(:) = 0._r8
! make lai average at landunit level
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       ttlai(p) = tlai(p)+tsai(p)
    enddo

    tlai_lu(:) = spval
    sumwt(:) = 0._r8
    do p = lbp,ubp
       if (ttlai(p) /= spval .and. wtlunit(p) /= 0._r8) then
          c = pcolumn(p)
          l = plandunit(p)
          if (sumwt(l) == 0._r8) tlai_lu(l) = 0._r8
          tlai_lu(l) = tlai_lu(l) + ttlai(p) * wtlunit(p)
          sumwt(l) = sumwt(l) + wtlunit(p)
       end if
    end do
    found = .false.
    do l = lbl,ubl
       if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = l
          exit
       else if (sumwt(l) /= 0._r8) then
          tlai_lu(l) = tlai_lu(l)/sumwt(l)
       end if
    end do
    if (found) then
       write(iulog,*) 'p2l_1d error: sumwt is greater than 1.0 at l= ',index
       call endrun()
    end if

! Loop through pfts

! initialize variables which get passed to the atmosphere
    flx_mss_vrt_dst(lbp:ubp,:)=0._r8

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       l = plandunit(p)
       
       ! the following code from subr. lnd_frc_mbl_get was adapted for lsm use
       ! purpose: return fraction of each gridcell suitable for dust mobilization
       
       ! the "bare ground" fraction of the current sub-gridscale cell decreases
       ! linearly from 1 to 0 as VAI(=tlai+tsai) increases from 0 to vai_mbl_thr
       ! if ice sheet, wetland, or lake, no dust allowed
       
       if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then
          if (tlai_lu(l) < vai_mbl_thr) then
             lnd_frc_mbl(p) = 1.0_r8 - (tlai_lu(l))/vai_mbl_thr
          else
             lnd_frc_mbl(p) = 0.0_r8
          endif
          lnd_frc_mbl(p) = lnd_frc_mbl(p) * (1.0_r8 - frac_sno(c))
       else          
          lnd_frc_mbl(p) = 0.0_r8   
       end if
    end do
    
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       if (lnd_frc_mbl(p)>1.0_r8 .or. lnd_frc_mbl(p)<0.0_r8) then
          write(iulog,*)'Error dstmbl: pft= ',p,' lnd_frc_mbl(p)= ',lnd_frc_mbl(p)
          call endrun 
       end if
    end do
    
    ! reset history output variables before next if-statement to avoid output = inf
    
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       flx_mss_vrt_dst_tot(p) = 0.0_r8
    end do
    do n = 1, ndst
       do fp = 1,num_nolakep
          p = filter_nolakep(fp)
          flx_mss_vrt_dst(p,n) = 0.0_r8
       end do
    end do
    
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       l = plandunit(p)
       g = pgridcell(p)
       
       ! only perform the following calculations if lnd_frc_mbl is non-zero 
       
       if (lnd_frc_mbl(p) > 0.0_r8) then

          ! the following comes from subr. frc_thr_rgh_fct_get
          ! purpose: compute factor by which surface roughness increases threshold
          !          friction velocity (currently a constant)
          
          frc_thr_rgh_fct = 1.0_r8
          
          ! the following comes from subr. frc_thr_wet_fct_get
          ! purpose: compute factor by which soil moisture increases threshold friction velocity
          ! adjust threshold velocity for inhibition by moisture
          ! modified 4/5/2002 (slevis) to use gravimetric instead of volumetric
          ! water content
          
          bd = (1._r8-watsat(c,1))*2.7e3_r8      ![kg m-3] Bulk density of dry surface soil
          gwc_sfc = h2osoi_vol(c,1)*SHR_CONST_RHOFW/bd    ![kg kg-1] Gravimetric H2O cont
          if (gwc_sfc > gwc_thr(c)) then
             frc_thr_wet_fct = sqrt(1.0_r8 + 1.21_r8 * (100.0_r8*(gwc_sfc - gwc_thr(c)))**0.68_r8)
          else
             frc_thr_wet_fct = 1.0_r8
          end if
          
          ! slevis: adding liqfrac here, because related to effects from soil water

          liqfrac = max( 0.0_r8, min( 1.0_r8, h2osoi_liq(c,1) / (h2osoi_ice(c,1)+h2osoi_liq(c,1)+1.0e-6_r8) ) )
          
          ! the following lines come from subr. dst_mbl
          ! purpose: adjust threshold friction velocity to acct for moisture and
          !          roughness. The ratio tmp1 / sqrt(forc_rho) comes from
          !          subr. wnd_frc_thr_slt_get which computes dry threshold
          !          friction velocity for saltation
          
          wnd_frc_thr_slt = tmp1 / sqrt(forc_rho(g)) * frc_thr_wet_fct * frc_thr_rgh_fct
          
          ! reset these variables which will be updated in the following if-block
          
          wnd_frc_slt = fv(p)
          flx_mss_hrz_slt_ttl = 0.0_r8
          flx_mss_vrt_dst_ttl(p) = 0.0_r8
          
          ! the following line comes from subr. dst_mbl
          ! purpose: threshold saltation wind speed
          
          wnd_rfr_thr_slt = u10(p) * wnd_frc_thr_slt / fv(p)
          
          ! the following if-block comes from subr. wnd_frc_slt_get 
          ! purpose: compute the saltating friction velocity
          ! theory: saltation roughens the boundary layer, AKA "Owen's effect"
          
          if (u10(p) >= wnd_rfr_thr_slt) then
             wnd_rfr_dlt = u10(p) - wnd_rfr_thr_slt
             wnd_frc_slt_dlt = 0.003_r8 * wnd_rfr_dlt * wnd_rfr_dlt
             wnd_frc_slt = fv(p) + wnd_frc_slt_dlt
          end if
          
          ! the following comes from subr. flx_mss_hrz_slt_ttl_Whi79_get
          ! purpose: compute vertically integrated streamwise mass flux of particles
          
          if (wnd_frc_slt > wnd_frc_thr_slt) then
             wnd_frc_rat = wnd_frc_thr_slt / wnd_frc_slt
             flx_mss_hrz_slt_ttl = cst_slt * forc_rho(g) * (wnd_frc_slt**3.0_r8) * &
                  (1.0_r8 - wnd_frc_rat) * (1.0_r8 + wnd_frc_rat) * (1.0_r8 + wnd_frc_rat) / grav
             
             ! the following loop originates from subr. dst_mbl
             ! purpose: apply land sfc and veg limitations and global tuning factor
             ! slevis: multiply flx_mss_hrz_slt_ttl by liqfrac to incude the effect 
             ! of frozen soil
             
             flx_mss_hrz_slt_ttl = flx_mss_hrz_slt_ttl * lnd_frc_mbl(p) * mbl_bsn_fct(c) * &
                  flx_mss_fdg_fct * liqfrac
          end if
          
          ! the following comes from subr. flx_mss_vrt_dst_ttl_MaB95_get
          ! purpose: diagnose total vertical mass flux of dust from vertically
          !          integrated streamwise mass flux
          
          dst_slt_flx_rat_ttl = 100.0_r8 * exp( log(10.0_r8) * (13.4_r8 * mss_frc_cly_vld(c) - 6.0_r8) )
          flx_mss_vrt_dst_ttl(p) = flx_mss_hrz_slt_ttl * dst_slt_flx_rat_ttl
          
       end if   ! lnd_frc_mbl > 0.0

    end do

    ! the following comes from subr. flx_mss_vrt_dst_prt in C. Zender's code
    ! purpose: partition total vertical mass flux of dust into transport bins
          
    do n = 1, ndst
       do m = 1, dst_src_nbr
          do fp = 1,num_nolakep
             p = filter_nolakep(fp)
             if (lnd_frc_mbl(p) > 0.0_r8) then
                flx_mss_vrt_dst(p,n) = flx_mss_vrt_dst(p,n) +  ovr_src_snk_mss(m,n) * flx_mss_vrt_dst_ttl(p)
             end if
          end do
       end do
    end do

    do n = 1, ndst
       do fp = 1,num_nolakep
          p = filter_nolakep(fp)
          if (lnd_frc_mbl(p) > 0.0_r8) then
             flx_mss_vrt_dst_tot(p) = flx_mss_vrt_dst_tot(p) + flx_mss_vrt_dst(p,n)
          end if
       end do
    end do

  end subroutine DustEmission

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine DustDryDep
!
! !INTERFACE:
!
  subroutine DustDryDep (lbp, ubp)
!
! !DESCRIPTION: 
!
! Determine Turbulent dry deposition for dust. Calculate the turbulent 
! component of dust dry deposition, (the turbulent deposition velocity 
! through the lowest atmospheric layer. CAM will calculate the settling 
! velocity through the whole atmospheric column. The two calculations 
! will determine the dust dry deposition flux to the surface.
! Note: Same process should occur over oceans. For the coupled CESM,
! we may find it more efficient to let CAM calculate the turbulent dep
! velocity over all surfaces. This would require passing the
! aerodynamic resistance, ram(1), and the friction velocity, fv, from
! the land to the atmosphere component. In that case, dustini need not
! calculate particle diamter (dmt_vwr) and particle density (dns_aer).
! Source: C. Zender's dry deposition code
!
! !USES
!
    use shr_const_mod, only : SHR_CONST_PI, SHR_CONST_RDAIR, SHR_CONST_BOLTZ
    use clm_atmlnd   , only : clm_a2l
!
! !ARGUMENTS:
!
    implicit none
    integer, intent(in) :: lbp, ubp     ! pft bounds
!
! !LOCAL VARIABLES
!
! local pointers to implicit in arguments
!
    integer , pointer :: plandunit(:)   ! pft's landunit index
    integer , pointer :: pgridcell(:)   ! pft's gridcell index
    integer , pointer :: ityplun(:)     ! landunit type
    real(r8), pointer :: pwtgcell(:)    ! weight of pft relative to corresponding gridcell
    real(r8), pointer :: forc_t(:)      ! atm temperature (K)
    real(r8), pointer :: forc_pbot(:)   ! atm pressure (Pa)
    real(r8), pointer :: forc_rho(:)    ! atm density (kg/m**3)
    real(r8), pointer :: fv(:)          ! friction velocity (m/s)
    real(r8), pointer :: ram1(:)        ! aerodynamical resistance (s/m)
    real(r8), pointer :: vlc_trb(:,:)   ! Turbulent deposn velocity (m/s)
    real(r8), pointer :: vlc_trb_1(:)   ! Turbulent deposition velocity 1
    real(r8), pointer :: vlc_trb_2(:)   ! Turbulent deposition velocity 2
    real(r8), pointer :: vlc_trb_3(:)   ! Turbulent deposition velocity 3
    real(r8), pointer :: vlc_trb_4(:)   ! Turbulent deposition velocity 4
!
! !REVISION HISTORY
! Created by Sam Levis
!
!
! !LOCAL VARIABLES
!EOP
!
    integer  :: p,l,g,m,n             ! indices
    real(r8) :: vsc_dyn_atm(lbp:ubp)  ! [kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(lbp:ubp)  ! [m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr_xpn           ! [frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
    real(r8) :: shm_nbr               ! [frc] Schmidt number
    real(r8) :: stk_nbr               ! [frc] Stokes number
    real(r8) :: mfp_atm               ! [m] Mean free path of air
    real(r8) :: dff_aer               ! [m2 s-1] Brownian diffusivity of particle
    real(r8) :: rss_trb               ! [s m-1] Resistance to turbulent deposition
    real(r8) :: slp_crc(lbp:ubp,ndst) ! [frc] Slip correction factor
    real(r8) :: vlc_grv(lbp:ubp,ndst) ! [m s-1] Settling velocity
    real(r8) :: rss_lmn(lbp:ubp,ndst) ! [s m-1] Quasi-laminar layer resistance
    real(r8) :: tmp                   ! temporary 

! constants

    real(r8),parameter::shm_nbr_xpn_lnd=-2._r8/3._r8 ![frc] shm_nbr_xpn over land
!------------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    forc_pbot => clm_a2l%forc_pbot
    forc_rho  => clm_a2l%forc_rho
    forc_t    => clm_a2l%forc_t

    ! Assign local pointers to derived type members (landunit-level)

    ityplun   => lun%itype

    ! Assign local pointers to derived type members (pft-level)

    plandunit => pft%landunit
    pgridcell => pft%gridcell
    pwtgcell  => pft%wtgcell  
    fv        => pps%fv
    ram1      => pps%ram1
    vlc_trb   => pdf%vlc_trb
    vlc_trb_1 => pdf%vlc_trb_1
    vlc_trb_2 => pdf%vlc_trb_2
    vlc_trb_3 => pdf%vlc_trb_3
    vlc_trb_4 => pdf%vlc_trb_4

    do p = lbp,ubp
       l = plandunit(p)
       ! Note: some glacier_mec pfts may have zero weight
       if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
          g = pgridcell(p)

          ! from subroutine dst_dps_dry (consider adding sanity checks from line 212)
          ! when code asks to use midlayer density, pressure, temperature,
          ! I use the data coming in from the atmosphere, ie forc_t, forc_pbot, forc_rho
          
          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          
          vsc_dyn_atm(p) = 1.72e-5_r8 * ((forc_t(g)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (forc_t(g)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm = 2.0_r8 * vsc_dyn_atm(p) / &   ![m] SeP97 p. 455
               (forc_pbot(g)*sqrt(8.0_r8/(SHR_CONST_PI*SHR_CONST_RDAIR*forc_t(g))))
          vsc_knm_atm(p) = vsc_dyn_atm(p) / forc_rho(g) ![m2 s-1] Kinematic viscosity of air
          
          do m = 1, ndst
             slp_crc(p,m) = 1.0_r8 + 2.0_r8 * mfp_atm * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
                  dmt_vwr(m)   ![frc] Slip correction factor SeP97 p. 464
             vlc_grv(p,m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
                  grav * slp_crc(p,m) / vsc_dyn_atm(p)   ![m s-1] Stokes' settling velocity SeP97 p. 466
             vlc_grv(p,m) = vlc_grv(p,m) * stk_crc(m)    ![m s-1] Correction to Stokes settling velocity
          end do
       end if
    end do

    do m = 1, ndst
       do p = lbp,ubp
          l = plandunit(p)
          if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
             g = pgridcell(p)
             
             stk_nbr = vlc_grv(p,m) * fv(p) * fv(p) / (grav * vsc_knm_atm(p))  ![frc] SeP97 p.965
             dff_aer = SHR_CONST_BOLTZ * forc_t(g) * slp_crc(p,m) / &          ![m2 s-1]
                  (3.0_r8*SHR_CONST_PI * vsc_dyn_atm(p) * dmt_vwr(m))          !SeP97 p.474
             shm_nbr = vsc_knm_atm(p) / dff_aer                                ![frc] SeP97 p.972
             shm_nbr_xpn = shm_nbr_xpn_lnd                                     ![frc]
             
             ! fxm: Turning this on dramatically reduces
             ! deposition velocity in low wind regimes
             ! Schmidt number exponent is -2/3 over solid surfaces and
             ! -1/2 over liquid surfaces SlS80 p. 1014
             ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd
             ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt # 
             
             tmp = shm_nbr**shm_nbr_xpn + 10.0_r8**(-3.0_r8/stk_nbr)
             rss_lmn(p,m) = 1.0_r8 / (tmp * fv(p)) ![s m-1] SeP97 p.972,965
          end if
       end do
    end do

    ! Lowest layer: Turbulent deposition (CAM will calc. gravitational dep)

    do m = 1, ndst
       do p = lbp,ubp
          l = plandunit(p)
          if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
             rss_trb = ram1(p) + rss_lmn(p,m) + ram1(p) * rss_lmn(p,m) * vlc_grv(p,m) ![s m-1]
             vlc_trb(p,m) = 1.0_r8 / rss_trb                                          ![m s-1]
          end if
       end do
    end do

    do p = lbp,ubp
       l = plandunit(p)
       if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
          vlc_trb_1(p) = vlc_trb(p,1)
          vlc_trb_2(p) = vlc_trb(p,2)
          vlc_trb_3(p) = vlc_trb(p,3)
          vlc_trb_4(p) = vlc_trb(p,4)
       end if
    end do 

  end subroutine DustDryDep

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine Dustini
!
! !INTERFACE:
!
  subroutine Dustini()
!
! !DESCRIPTION: 
!
! Compute source efficiency factor from topography
! Initialize other variables used in subroutine Dust:
! ovr_src_snk_mss(m,n) and tmp1.
! Define particle diameter and density needed by atm model
! as well as by dry dep model
! Source: Paul Ginoux (for source efficiency factor)
! Modifications by C. Zender and later by S. Levis
! Rest of subroutine from C. Zender's dust model
!
! !USES
    use shr_const_mod, only: SHR_CONST_PI, SHR_CONST_RDAIR
    use decompMod, only : get_proc_bounds
    use shr_spfn_mod, only: erf => shr_spfn_erf
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY
! Created by Samual Levis
!
! !LOCAL VARIABLES
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: mbl_bsn_fct(:) !basin factor
!
!
! !LOCAL VARIABLES
!EOP
!
    integer  :: fc,c,l,m,n              ! indices
    real(r8) :: ovr_src_snk_frc
    real(r8) :: sqrt2lngsdi             ! [frc] Factor in erf argument
    real(r8) :: lndmaxjovrdmdni         ! [frc] Factor in erf argument
    real(r8) :: lndminjovrdmdni         ! [frc] Factor in erf argument
    real(r8) :: ryn_nbr_frc_thr_prx_opt ! [frc] Threshold friction Reynolds number approximation for optimal size
    real(r8) :: ryn_nbr_frc_thr_opt_fnc ! [frc] Threshold friction Reynolds factor for saltation calculation
    real(r8) :: icf_fct                 ! Interpartical cohesive forces factor for saltation calc
    real(r8) :: dns_fct                 ! Density ratio factor for saltation calculation
    real(r8) :: dmt_min(ndst)           ! [m] Size grid minimum
    real(r8) :: dmt_max(ndst)           ! [m] Size grid maximum
    real(r8) :: dmt_ctr(ndst)           ! [m] Diameter at bin center
    real(r8) :: dmt_dlt(ndst)           ! [m] Width of size bin
    real(r8) :: slp_crc(ndst)           ! [frc] Slip correction factor
    real(r8) :: vlm_rsl(ndst)           ! [m3 m-3] Volume concentration resolved
    real(r8) :: vlc_stk(ndst)           ! [m s-1] Stokes settling velocity
    real(r8) :: vlc_grv(ndst)           ! [m s-1] Settling velocity
    real(r8) :: ryn_nbr_grv(ndst)       ! [frc] Reynolds number at terminal velocity
    real(r8) :: cff_drg_grv(ndst)       ! [frc] Drag coefficient at terminal velocity
    real(r8) :: tmp                     ! temporary 
    real(r8) :: ln_gsd                  ! [frc] ln(gsd)
    real(r8) :: gsd_anl                 ! [frc] Geometric standard deviation
    real(r8) :: dmt_vma                 ! [m] Mass median diameter analytic She84 p.75 Tabl.1
    real(r8) :: dmt_nma                 ! [m] Number median particle diameter
    real(r8) :: lgn_dst                 ! Lognormal distribution at sz_ctr
    real(r8) :: eps_max                 ! [frc] Relative accuracy for convergence
    real(r8) :: eps_crr                 ! [frc] Current relative accuracy
    real(r8) :: itr_idx                 ! [idx] Counting index
    real(r8) :: dns_mdp                 ! [kg m-3] Midlayer density
    real(r8) :: mfp_atm                 ! [m] Mean free path of air
    real(r8) :: vsc_dyn_atm             ! [kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm             ! [kg m-1 s-1] Kinematic viscosity of air
    real(r8) :: vlc_grv_old             ! [m s-1] Previous gravitational settling velocity
    real(r8) :: series_ratio            ! Factor for logarithmic grid
    real(r8) :: lngsdsqrttwopi_rcp      ! Factor in lognormal distribution
    real(r8) :: sz_min(sz_nbr)          ! [m] Size Bin minima
    real(r8) :: sz_max(sz_nbr)          ! [m] Size Bin maxima
    real(r8) :: sz_ctr(sz_nbr)          ! [m] Size Bin centers
    real(r8) :: sz_dlt(sz_nbr)          ! [m] Size Bin widths
    
    ! constants
    real(r8) :: dmt_vma_src(dst_src_nbr) =    &     ! [m] Mass median diameter
         (/ 0.832e-6_r8 , 4.82e-6_r8 , 19.38e-6_r8 /)        ! BSM96 p. 73 Table 2
    real(r8) :: gsd_anl_src(dst_src_nbr) =    &     ! [frc] Geometric std deviation
         (/ 2.10_r8     ,  1.90_r8   , 1.60_r8     /)        ! BSM96 p. 73 Table 2
    real(r8) :: mss_frc_src(dst_src_nbr) =    &     ! [frc] Mass fraction 
         (/ 0.036_r8, 0.957_r8, 0.007_r8 /)                  ! BSM96 p. 73 Table 2
    real(r8) :: dmt_grd(5) =                  &     ! [m] Particle diameter grid
         (/ 0.1e-6_r8, 1.0e-6_r8, 2.5e-6_r8, 5.0e-6_r8, 10.0e-6_r8 /)
    real(r8), parameter :: dmt_slt_opt = 75.0e-6_r8    ! [m] Optim diam for saltation
    real(r8), parameter :: dns_slt = 2650.0_r8         ! [kg m-3] Density of optimal saltation particles

    ! declare erf intrinsic function
    real(r8) :: dum     ! dummy variable for erf test

    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices 
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
!------------------------------------------------------------------------

    ! Assign local pointers to derived type scalar members (column-level)

    mbl_bsn_fct => cps%mbl_bsn_fct

    ! the following comes from (1) szdstlgn.F subroutine ovr_src_snk_frc_get
    !                      and (2) dstszdst.F subroutine dst_szdst_ini
    ! purpose(1): given one set (the "source") of lognormal distributions,
    !             and one set of bin boundaries (the "sink"), compute and return
    !             the overlap factors between the source and sink distributions
    ! purpose(2): set important statistics of size distributions

    do m = 1, dst_src_nbr
       sqrt2lngsdi = sqrt(2.0_r8) * log(gsd_anl_src(m))
       do n = 1, ndst
          lndmaxjovrdmdni = log(dmt_grd(n+1)/dmt_vma_src(m))
          lndminjovrdmdni = log(dmt_grd(n  )/dmt_vma_src(m))
          ovr_src_snk_frc = 0.5_r8 * (erf(lndmaxjovrdmdni/sqrt2lngsdi) - &
                                   erf(lndminjovrdmdni/sqrt2lngsdi))
          ovr_src_snk_mss(m,n) = ovr_src_snk_frc * mss_frc_src(m)
       end do
    end do

    ! The following code from subroutine wnd_frc_thr_slt_get was placed 
    ! here because tmp1 needs to be defined just once

    ryn_nbr_frc_thr_prx_opt = 0.38_r8 + 1331.0_r8 * (100.0_r8*dmt_slt_opt)**1.56_r8

    if (ryn_nbr_frc_thr_prx_opt < 0.03_r8) then
       write(iulog,*) 'dstmbl: ryn_nbr_frc_thr_prx_opt < 0.03'
       call endrun
    else if (ryn_nbr_frc_thr_prx_opt < 10.0_r8) then
       ryn_nbr_frc_thr_opt_fnc = -1.0_r8 + 1.928_r8 * (ryn_nbr_frc_thr_prx_opt**0.0922_r8)
       ryn_nbr_frc_thr_opt_fnc = 0.1291_r8 * 0.1291_r8 / ryn_nbr_frc_thr_opt_fnc
    else
       ryn_nbr_frc_thr_opt_fnc = 1.0_r8 - 0.0858_r8 * exp(-0.0617_r8*(ryn_nbr_frc_thr_prx_opt-10.0_r8))
       ryn_nbr_frc_thr_opt_fnc = 0.120_r8 * 0.120_r8 * ryn_nbr_frc_thr_opt_fnc * ryn_nbr_frc_thr_opt_fnc
    end if

    icf_fct = 1.0_r8 + 6.0e-07_r8 / (dns_slt * grav * (dmt_slt_opt**2.5_r8))
    dns_fct = dns_slt * grav * dmt_slt_opt
    tmp1 = sqrt(icf_fct * dns_fct * ryn_nbr_frc_thr_opt_fnc)

    ! Set basin factor to 1 for now

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    do c = begc, endc
      l = col%landunit(c)
      if (.not. lun%lakpoi(l)) then
         mbl_bsn_fct(c) = 1.0_r8
      end if
    end do

    ! Introducing particle diameter. Needed by atm model and by dry dep model.
    ! Taken from Charlie Zender's subroutines dst_psd_ini, dst_sz_rsl,
    ! grd_mk (dstpsd.F90) and subroutine lgn_evl (psdlgn.F90)
    
    ! Charlie allows logarithmic or linear option for size distribution
    ! however, he hardwires the distribution to logarithmic in his code
    ! therefore, I take his logarithmic code only
    ! furthermore, if dst_nbr == 4, he overrides the automatic grid calculation
    ! he currently works with dst_nbr = 4, so I only take the relevant code
    ! if ndst ever becomes different from 4, must add call grd_mk (dstpsd.F90)
    ! as done in subroutine dst_psd_ini
    ! note that here ndst = dst_nbr
    
    ! Override automatic grid with preset grid if available

    if (ndst == 4) then
       do n = 1, ndst
          dmt_min(n) = dmt_grd(n)                       ![m] Max diameter in bin
          dmt_max(n) = dmt_grd(n+1)                     ![m] Min diameter in bin
          dmt_ctr(n) = 0.5_r8 * (dmt_min(n)+dmt_max(n)) ![m] Diameter at bin ctr
          dmt_dlt(n) = dmt_max(n)-dmt_min(n)            ![m] Width of size bin
       end do
    else
       write(iulog,*) 'Dustini error: ndst must equal to 4 with current code'
       call endrun                                      !see more comments above
    end if                                              !end if ndst == 4

    ! Bin physical properties

    gsd_anl = 2.0_r8      ! [frc] Geometric std dev PaG77 p. 2080 Table1
    ln_gsd = log(gsd_anl)
    dns_aer = 2.5e+3_r8   ! [kg m-3] Aerosol density

    ! Set a fundamental statistic for each bin

    dmt_vma = 3.5000e-6_r8 ! [m] Mass median diameter analytic She84 p.75 Table1

    ! Compute analytic size statistics
    ! Convert mass median diameter to number median diameter (call vma2nma)

    dmt_nma = dmt_vma * exp(-3.0_r8*ln_gsd*ln_gsd) ! [m]

    ! Compute resolved size statistics for each size distribution
    ! In C. Zender's code call dst_sz_rsl

    do n = 1, ndst

       series_ratio = (dmt_max(n)/dmt_min(n))**(1.0_r8/sz_nbr)
       sz_min(1) = dmt_min(n)
       do m = 2, sz_nbr                            ! Loop starts at 2
          sz_min(m) = sz_min(m-1) * series_ratio
       end do

       ! Derived grid values
       do m = 1, sz_nbr-1                          ! Loop ends at sz_nbr-1
          sz_max(m) = sz_min(m+1)                  ! [m]
       end do
       sz_max(sz_nbr) = dmt_max(n)                 ! [m]

       ! Final derived grid values
       do m = 1, sz_nbr
          sz_ctr(m) = 0.5_r8 * (sz_min(m)+sz_max(m))
          sz_dlt(m) = sz_max(m)-sz_min(m)
       end do

       lngsdsqrttwopi_rcp = 1.0_r8 / (ln_gsd*sqrt(2.0_r8*SHR_CONST_PI))
       dmt_vwr(n) = 0.0_r8 ! [m] Mass wgted diameter resolved
       vlm_rsl(n) = 0.0_r8 ! [m3 m-3] Volume concentration resolved

       do m = 1, sz_nbr

          ! Evaluate lognormal distribution for these sizes (call lgn_evl)
          tmp = log(sz_ctr(m)/dmt_nma) / ln_gsd
          lgn_dst = lngsdsqrttwopi_rcp * exp(-0.5_r8*tmp*tmp) / sz_ctr(m)

          ! Integrate moments of size distribution
          dmt_vwr(n) = dmt_vwr(n) + sz_ctr(m) *                    &
               SHR_CONST_PI / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
          vlm_rsl(n) = vlm_rsl(n) +                                &
               SHR_CONST_PI / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn

       end do

       dmt_vwr(n) = dmt_vwr(n) / vlm_rsl(n) ![m] Mass weighted diameter resolved

    end do

    ! calculate correction to Stokes' settling velocity (subroutine stk_crc_get)

    eps_max = 1.0e-4_r8
    dns_mdp = 100000._r8 / (295.0_r8*SHR_CONST_RDAIR) ![kg m-3] const prs_mdp & tpt_vrt

    ! Size-independent thermokinetic properties

    vsc_dyn_atm = 1.72e-5_r8 * ((295.0_r8/273.0_r8)**1.5_r8) * 393.0_r8 / &
         (295.0_r8+120.0_r8)      ![kg m-1 s-1] RoY94 p.102 tpt_mdp=295.0
    mfp_atm = 2.0_r8 * vsc_dyn_atm / &  !SeP97 p. 455 constant prs_mdp, tpt_mdp
         (100000._r8*sqrt(8.0_r8/(SHR_CONST_PI*SHR_CONST_RDAIR*295.0_r8)))
    vsc_knm_atm = vsc_dyn_atm / dns_mdp ![m2 s-1] Kinematic viscosity of air

    do m = 1, ndst
       slp_crc(m) = 1.0_r8 + 2.0_r8 * mfp_atm *                      &
            (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
            dmt_vwr(m)                      ! [frc] Slip correction factor SeP97 p.464
       vlc_stk(m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
            grav * slp_crc(m) / vsc_dyn_atm ! [m s-1] SeP97 p.466
    end do

    ! For Reynolds number flows Re < 0.1 Stokes' velocity is valid for
    ! vlc_grv SeP97 p. 466 (8.42). For larger Re, inertial effects become
    ! important and empirical drag coefficients must be employed
    ! Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
    ! Using Stokes' velocity rather than iterative solution with empirical
    ! drag coefficient causes 60% errors for D = 200 um SeP97 p. 468

    ! Iterative solution for drag coefficient, Reynolds number, and terminal veloc
    do m = 1, ndst

       ! Initialize accuracy and counter
       eps_crr = eps_max + 1.0_r8  ![frc] Current relative accuracy
       itr_idx = 0                 ![idx] Counting index

       ! Initial guess for vlc_grv is exact for Re < 0.1
       vlc_grv(m) = vlc_stk(m)     ![m s-1]

       do while(eps_crr > eps_max)

          ! Save terminal velocity for convergence test
          vlc_grv_old = vlc_grv(m) ![m s-1]
          ryn_nbr_grv(m) = vlc_grv(m) * dmt_vwr(m) / vsc_knm_atm !SeP97 p.460

          ! Update drag coefficient based on new Reynolds number
          if (ryn_nbr_grv(m) < 0.1_r8) then
             cff_drg_grv(m) = 24.0_r8 / ryn_nbr_grv(m) !Stokes' law Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0_r8) then
             cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) *    &
                  (1.0_r8 + 3.0_r8*ryn_nbr_grv(m)/16.0_r8 + &
                  9.0_r8*ryn_nbr_grv(m)*ryn_nbr_grv(m)*     &
                  log(2.0_r8*ryn_nbr_grv(m))/160.0_r8)        !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 500.0_r8) then
             cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) * &
                  (1.0_r8 + 0.15_r8*ryn_nbr_grv(m)**0.687_r8) !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0e5_r8) then
             cff_drg_grv(m) = 0.44_r8                         !Sep97 p.463 (8.32)
          else
             write(iulog,'(a,es9.2)') "ryn_nbr_grv(m) = ",ryn_nbr_grv(m)
             write(iulog,*)'Dustini error: Reynolds number too large in stk_crc_get()'
             call endrun 
          end if

          ! Update terminal velocity based on new Reynolds number and drag coeff
          ! [m s-1] Terminal veloc SeP97 p.467 (8.44)

          vlc_grv(m) = sqrt(4.0_r8 * grav * dmt_vwr(m) * slp_crc(m) * dns_aer / &
               (3.0_r8*cff_drg_grv(m)*dns_mdp))   
          eps_crr = abs((vlc_grv(m)-vlc_grv_old)/vlc_grv(m)) !Relative convergence
          if (itr_idx == 12) then
             ! Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0
             ! due to discontinuities in derivative of drag coefficient
             vlc_grv(m) = 0.5_r8 * (vlc_grv(m)+vlc_grv_old)  ! [m s-1]
          end if
          if (itr_idx > 20) then
             write(iulog,*) 'Dustini error: Terminal velocity not converging ',&
                  ' in stk_crc_get(), breaking loop...'
             goto 100                                        !to next iteration
          end if
          itr_idx = itr_idx + 1

       end do                                                !end while

100    continue   !Label to jump to when iteration does not converge
    end do   !end loop over size

    ! Compute factors to convert Stokes' settling velocities to
    ! actual settling velocities

    do m = 1, ndst
       stk_crc(m) = vlc_grv(m) / vlc_stk(m)
    end do

  end subroutine Dustini

end module DUSTMod
