
      module mo_strato_rates
!=======================================================================
! ROUTINE
!   ratecon_sfstrat.f
!
!  Date...
!  15 August 2002
!  11 April  2008
!
!  Programmed by...
!   Douglas E. Kinnison
!
! DESCRIPTION
!
! Derivation of the rate constant for reactions on
!   sulfate, NAT, and ICE aerosols.
!
!
! Sulfate Aerosol Reactions               Rxn#   Gamma
!   N2O5   + H2O(l)     =>  2HNO3         (1)    f(wt%)
!   ClONO2 + H2O(l)     =>  HOCl + HNO3   (2)    f(T,P,HCl,H2O,r)
!   BrONO2 + H2O(l)     =>  HOBr + HNO3   (3)    f(T,P,H2O,r)
!   ClONO2 + HCl(l)     =>  Cl2  + HNO3   (4)    f(T,P,HCl,H2O,r)
!   HOCl   + HCl(l)     =>  Cl2  + H2O    (5)    f(T,P,HCl,HOCl,H2O,r)
!   HOBr   + HCl(l)     =>  BrCl + H2O    (6)    f(T,P,HCl,HOBr,H2O,r)
!
! Nitric Acid Di-hydrate Reactions        Rxn#    Gamma   Reference
!   N2O5   + H2O(s)     =>  2HNO3         (7)     4e-4   JPL06-2
!   ClONO2 + H2O(s)     =>  HOCl + HNO3   (8)     4e-3   JPL06-2
!   ClONO2 + HCl(s)     =>  Cl2  + HNO3   (9)     0.2    JPL06-2
!   HOCl   + HCl(s)     =>  Cl2  + H2O    (10)    0.1    JPL06-2
!   BrONO2 + H2O(s)     =>  HOBr + HNO3   (11)    0.3    David Hanson PC
!
! ICE Aersol Reactions                    Rxn#    Gamma
!   N2O5   + H2O(s)     =>  2HNO3         (12)     0.02   JPL06-2
!   ClONO2 + H2O(s)     =>  HOCl + HNO3   (13)     0.3    JPL06-2
!   BrONO2 + H2O(s)     =>  HOBr + HNO3   (14)     0.3    JPL06-2
!   ClONO2 + HCl(s)     =>  Cl2  + HNO3   (15)     0.3    JPL06-2
!   HOCl   + HCl(s)     =>  Cl2  + H2O    (16)     0.2    JPL06-2
!   HOBr   + HCl(s)     =>  BrCl + H2O    (17)     0.3    JPL06-2
!
! NOTE: The rate constants derived from species reacting with H2O are
!       first order (i.e., sec-1 units) - an example is N2O5 + H2O = 2HNO3.
!       Other reactions, e.g., ClONO2 + HCl have rate constants that
!       are second order (i.e., cm+3 molecules-1 sec-1 units). In all
!       of these types of reactions the derived first order rate constant
!       {0.25*(mean Velocity)*SAD*gamma} is divided by the HCl abundance
!       to derive the correct second order units.
!
! NOTE: Liquid Sulfate Aerosols...
!       See coding for references on how the Sulfate Aerosols were handled.
!       Data was used that was more recent than JPL00.
!
!
! INPUT:
!  ad      .    .... air density, molec. cm-3
!  pmid        ..... pressures, hPa
!  temp        ..... temperatures, K
!  rad_sulfate ..... Surface area density, cm2 cm-3
!  sad_sulfate ..... Surface area density, cm2 cm-3
!  sad_nat     ..... Surface area density, cm2 cm-3
!  sad_ice     ..... Surface area density, cm2 cm-3
!  brono2mv    ..... BrONO2 Volume Mixing Ratio
!  clono2mvr   ..... ClONO2 Volume Mixing Ratio
!  h2omvr      ..... H2O Volume Mixing Ratio
!  hclmvr      ..... HCl Volume Mixing Ratio
!  hobrmvr     ..... HOBr Volume Mixing Ratio
!  hoclmvr     ..... HOCl Volume Mixing Ratio
!  n2o5mvr     ..... N2O5 Volume Mixing Ratio
!
! OUTPUT:
!
!  rxt         ..... Rate constant (s-1 and cm3 sec-1 molec-1)
!=======================================================================

      private
      public :: ratecon_sfstrat, init_strato_rates, has_strato_chem

      integer :: id_brono2, id_clono2, id_hcl, id_hocl, &
           id_hobr, id_n2o5
      integer :: rid_het1,  rid_het2,  rid_het3,  rid_het4,  rid_het5, &
           rid_het6,  rid_het7,  rid_het8,  rid_het9,  rid_het10, &
           rid_het11, rid_het12, rid_het13, rid_het14, rid_het15, &
           rid_het16, rid_het17

      logical :: has_strato_chem 

      contains

        subroutine init_strato_rates

          use mo_chem_utls, only : get_rxt_ndx, get_spc_ndx
          use mo_aero_settling, only: strat_aer_settl_init
          implicit none

          integer :: ids(23)

          rid_het1  = get_rxt_ndx( 'het1' )
          rid_het2  = get_rxt_ndx( 'het2' )
          rid_het3  = get_rxt_ndx( 'het3' )
          rid_het4  = get_rxt_ndx( 'het4' )
          rid_het5  = get_rxt_ndx( 'het5' )
          rid_het6  = get_rxt_ndx( 'het6' )
          rid_het7  = get_rxt_ndx( 'het7' )
          rid_het8  = get_rxt_ndx( 'het8' )
          rid_het9  = get_rxt_ndx( 'het9' )
          rid_het10 = get_rxt_ndx( 'het10' )
          rid_het11 = get_rxt_ndx( 'het11' )
          rid_het12 = get_rxt_ndx( 'het12' )
          rid_het13 = get_rxt_ndx( 'het13' )
          rid_het14 = get_rxt_ndx( 'het14' )
          rid_het15 = get_rxt_ndx( 'het15' )
          rid_het16 = get_rxt_ndx( 'het16' )
          rid_het17 = get_rxt_ndx( 'het17' )

          id_brono2 = get_spc_ndx( 'BRONO2' )
          id_clono2 = get_spc_ndx( 'CLONO2' )
          id_hcl    = get_spc_ndx( 'HCL' )
          id_hocl   = get_spc_ndx( 'HOCL' )
          id_hobr   = get_spc_ndx( 'HOBR' )
          id_n2o5   = get_spc_ndx( 'N2O5' )

          ids(:) = (/ rid_het1, rid_het2, rid_het3, rid_het4, rid_het5, rid_het6, rid_het7, rid_het8, &
               rid_het9, rid_het10, rid_het11, rid_het12, rid_het13, rid_het14, rid_het15, &
               rid_het16, rid_het17, id_brono2, id_clono2, id_hcl, id_hocl, id_hobr, id_n2o5 /)

          has_strato_chem = all( ids(:) > 0 )

          if (.not. has_strato_chem) return

          call strat_aer_settl_init

        endsubroutine init_strato_rates

      subroutine ratecon_sfstrat( ad, pmid, temp, rad_sulfate, sad_sulfate, &
                                  sad_nat, sad_ice, h2ovmr, vmr, rxt, ncol )

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods,    only : adv_mass, rxntot, gas_pcnst
      use ppgrid,       only : pcols, pver
      use mo_sad,       only : sad_top
      use cam_logfile,  only : iulog

      implicit none

!-----------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: ncol                                 ! columns in chunk
      real(r8), dimension(ncol,pver,gas_pcnst), intent(in) :: &   ! species concentrations (mol/mol)
        vmr
      real(r8), dimension(ncol,pver), intent(in) :: &
        ad, &                                                     ! Air Density (molec. cm-3)
        rad_sulfate, &                                            ! Radius of Sulfate Aerosol (cm)
        sad_ice, &                                                ! ICE Surface Area Density (cm-1)
        sad_nat, &                                                ! NAT Surface Area Density (cm-1)
        sad_sulfate, &                                            ! Sulfate Surface Area Density (cm-1)
        h2ovmr                                                    ! water vapor volume mixing ratio( gas phase )
      real(r8), dimension(pcols,pver), intent(in) :: &
        pmid, &                                                   ! pressure (Pa)
        temp                                                      ! temperature (K)

      real(r8), intent(out) :: &
        rxt(ncol,pver,rxntot)                                     ! rate constants

!-----------------------------------------------------------------------
!  	... local variables
!-----------------------------------------------------------------------
	real(r8), parameter :: small_conc = 1.e-30_r8
	real(r8), parameter :: av_const   = 2.117265e4_r8  ! (8*8.31448*1000 / PI)
	real(r8), parameter :: pa2mb      = 1.e-2_r8       ! Pa to mb
	real(r8), parameter :: m2cm       = 100._r8        ! meters to cms

      integer :: &
        i, &                      ! altitude loop index
        k, &                      ! level loop index
        m                         ! species index

!-----------------------------------------------------------------------
!   	... variables for gamma calculations
!-----------------------------------------------------------------------
      real(r8) :: &
        brono2vmr, &                            ! BrONO2 Volume Mixing Ratio
        clono2vmr, &                            ! ClONO2 Volume Mixing Ratio
        hclvmr, &                               ! HCl Volume Mixing Ratio
        hcldeni, &                              ! inverse of HCl density
        cntdeni, &                              ! inverse of ClONO2 density
        hocldeni, &                             ! inverse of HOCl density
        hobrdeni, &                             ! inverse of HOBr density
        hoclvmr, &                              ! HOCl Volume Mixing Ratio
        hobrvmr, &                              ! HOBr Volume Mixing Ratio
        n2o5vmr                                 ! N2O5 Volume Mixing Ratio

        real(r8) :: &
        av_n2o5, &                              ! N2O5 Mean Velocity (cm s-1)
        av_clono2, &                            ! ClONO2 Mean Velocity (cm s-1)
        av_brono2, &                            ! BrONO2Mean Velocity (cm s-1)
        av_hocl, &                              ! HOCl Mean Velocity (cm s-1)
        av_hobr                                 ! HOBr Mean Velocity (cm s-1)

      real(r8) :: &
        pzero_h2o, &                            ! H2O sat vapor press (mbar)
        e0, e1, e2, e3, &                       ! coefficients for H2O sat vapor press.
        aw, &                                   ! Water activity
        m_h2so4, &                              ! H2SO4 molality (mol/kg)
        wt, &                                   ! wt % H2SO4
        y1, y2, &                               ! used in H2SO4 molality
        &  a1, b1, c1, d1, &                    ! used in H2SO4 molality
        a2, b2, c2, d2                          ! used in H2SO4 molality

        real(r8) :: &
        z1, z2, z3, &                           ! used in H2SO4 soln density
        den_h2so4, &                            ! H2SO4 soln density, g/cm3
        mol_h2so4, &                            ! Molality of H2SO4, mol / kg
        molar_h2so4, &                          ! Molarity of H2SO4, mol / l
        x_h2so4, &                              ! H2SO4 mole fraction
        aconst, tzero, &                        ! used in viscosity of H2SO4
        vis_h2so4, &                            ! H2SO4 viscosity
        ah, &                                   ! Acid activity, molarity units
        term1,term2,term3,term4, &              ! used in ah
        term5,term6,term7,term0, &
        T_limit, &                              ! temporary variable for temp (185-260K range)
        T_limiti, &                             ! 1./T_limit
        T_limitsq, &                            ! sqrt( T_limit )
        rad_sulf, &                             ! temporary variable for sulfate radius (cm)
        sadsulf, &                              ! temporary variable for sulfate radius (cm)
        sadice, &                               ! temporary variable for sulfate radius (cm)
        sadnat                                  ! temporary variable for sulfate radius (cm)

      real(r8) :: &
        C_cnt, S_cnt, &                         ! used in H_cnt
        H_cnt, &                                ! Henry's law coeff. for ClONO2
        H_hcl, &                                ! Henry's law coeff. for HCl
        D_cnt, &
        k_hydr, &
        k_h2o, &
        k_h, &
        k_hcl, &
        rdl_cnt, &
        f_cnt, &
        M_hcl, &
        atmos

      real(r8) :: &
        Gamma_b_h2o, &
        Gamma_cnt_rxn, &
        Gamma_b_hcl, &
        Gamma_s, &
        Fhcl, &
        Gamma_s_prime, &
        Gamma_b_hcl_prime, &
        Gamma_b, &
        gprob_n2o5, &
        gprob_rxn, &
        gprob_tot, &
        gprob_cnt, &
        gprob_cnt_hcl, &
        gprob_cnt_h2o

        real(r8) :: &
        D_hocl, &
        k_hocl_hcl, &
        C_hocl, &
        S_hocl, &
        H_hocl, &
        Gamma_hocl_rxn, &
        rdl_hocl, &
        f_hocl, &
        gprob_hocl_hcl

        real(r8) :: &
        h1, h2, h3, &
        alpha, &
        gprob_bnt_h2o

      real(r8) :: &
        C_hobr, &
        D_hobr, &
        aa, bb, cc, dd, &
        k_hobr_hcl, &
        k_dl, &
        k_wasch, &
        H_hobr, &
        rdl_hobr, &
        Gamma_hobr_rxn, &
        f_hobr, &
        gprob_hobr_hcl

      real(r8) :: &
        pmb,&					! Pressure, mbar (hPa)
        pH2O_atm,&				! Partial press. H2O (atm)
        pH2O_hPa,&				! Partial press. H2O (hPa)
        pHCl_atm,&				! Partial press. HCl (atm)
        pCNT_atm                                ! Partial press. ClONO2 (atm)

!-----------------------------------------------------------------------
!     	... Used in pzero h2o calculation
!-----------------------------------------------------------------------
      real(r8), parameter :: wt_e0 = 18.452406985_r8
      real(r8), parameter :: wt_e1 = 3505.1578807_r8
      real(r8), parameter :: wt_e2 = 330918.55082_r8
      real(r8), parameter :: wt_e3 = 12725068.262_r8

      real(r8) :: &
        wrk, tmp

      real(r8), parameter :: small = 1.e-16_r8

      if (.not. has_strato_chem) return

!-----------------------------------------------------------------------
!     	... intialize rate constants
!-----------------------------------------------------------------------
      do k = 1,pver
         rxt(:,k,rid_het1) = 0._r8
         rxt(:,k,rid_het2) = 0._r8
         rxt(:,k,rid_het3) = 0._r8
         rxt(:,k,rid_het4) = 0._r8
         rxt(:,k,rid_het5) = 0._r8
         rxt(:,k,rid_het6) = 0._r8
         rxt(:,k,rid_het7) = 0._r8
         rxt(:,k,rid_het8) = 0._r8
         rxt(:,k,rid_het9) = 0._r8
         rxt(:,k,rid_het10) = 0._r8
         rxt(:,k,rid_het11) = 0._r8
         rxt(:,k,rid_het12) = 0._r8
         rxt(:,k,rid_het13) = 0._r8
         rxt(:,k,rid_het14) = 0._r8
         rxt(:,k,rid_het15) = 0._r8
         rxt(:,k,rid_het16) = 0._r8
         rxt(:,k,rid_het17) = 0._r8
      end do

!-----------------------------------------------------------------------
!     	... set rate constants
!-----------------------------------------------------------------------
Level_loop : &
      do k = sad_top+1,pver
column_loop : &
         do i = 1,ncol
!-----------------------------------------------------------------------
!	... set species, pmb, and atmos
!-----------------------------------------------------------------------
	    brono2vmr = vmr(i,k,id_brono2)
	    clono2vmr = vmr(i,k,id_clono2)
	    hclvmr    = vmr(i,k,id_hcl)
	    hoclvmr      = vmr(i,k,id_hocl)
	    hobrvmr      = vmr(i,k,id_hobr)
	    if( hclvmr > 0._r8 ) then
	       hcldeni  = 1._r8/(hclvmr*ad(i,k))
	    end if
	    if( clono2vmr > 0._r8 ) then
	       cntdeni  = 1._r8/(clono2vmr*ad(i,k))
	    end if
	    if( hoclvmr > 0._r8 ) then
	       hocldeni  = 1._r8/(hoclvmr*ad(i,k))
	    end if
	    if( hobrvmr > 0._r8 ) then
	       hobrdeni  = 1._r8/(hobrvmr*ad(i,k))
	    end if
	    n2o5vmr      = vmr(i,k,id_n2o5)
            sadsulf      = sad_sulfate(i,k)
            sadnat       = sad_nat(i,k)
	    sadice       = sad_ice(i,k)
            pmb          = pa2mb*pmid(i,k)
            atmos        = pmb/1013.25_r8

!-----------------------------------------------------------------------
!  	... setup for stratospheric aerosols
!           data range set: 185K - 240K;    GRL, 24, 1931, 1997
!-----------------------------------------------------------------------
            T_limit   = max( temp(i,k),185._r8 )
            T_limit   = min( T_limit,240._r8 )
            T_limiti  = 1._r8/T_limit
            T_limitsq = sqrt( T_limit )

!-----------------------------------------------------------------------
!     .... Average velocity (8RT*1000/(PI*MW))**1/2 * 100.(units cm s-1)
!     .... or (av_const*T/M2)**1/2
!-----------------------------------------------------------------------
	    wrk       = av_const*T_limit
            av_n2o5   = sqrt( wrk/adv_mass(id_n2o5) )*m2cm
            av_clono2 = sqrt( wrk/adv_mass(id_clono2) )*m2cm
            av_brono2 = sqrt( wrk/adv_mass(id_brono2) )*m2cm
            av_hocl   = sqrt( wrk/adv_mass(id_hocl) )*m2cm
            av_hobr   = sqrt( wrk/adv_mass(id_hobr) )*m2cm
has_sadsulf : &
            if( sadsulf > 0._r8 ) then
!-----------------------------------------------------------------------
!     .... Partial Pressure of H2O, ClONO2, and HCl in atmospheres
!-----------------------------------------------------------------------
               if( hclvmr > 0._r8 ) then
                  pHCl_atm  = hclvmr*atmos
               else
                  pHCl_atm  = 0._r8
               end if

               if( clono2vmr > 0._r8 ) then
                  pCNT_atm  = clono2vmr*atmos
               else
                  pCNT_atm  = 0._r8
               end if

               if( h2ovmr(i,k) > 0._r8 ) then
                  pH2O_atm  = h2ovmr(i,k)*atmos
               else
                  pH2O_atm  = 0._r8
               end if

!-----------------------------------------------------------------------
!     .... Partial Pressure of H2O in hPa
!-----------------------------------------------------------------------
               pH2O_hpa = h2ovmr(i,k)*pmb
!-----------------------------------------------------------------------
!     .... Calculate the h2so4 Wt% and Activity of H2O - 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     ... Saturation Water Vapor Pressure (mbar)
!-----------------------------------------------------------------------
               pzero_h2o = exp( wt_e0 - T_limiti*(wt_e1 + T_limiti*(wt_e2 - T_limiti*wt_e3)) )
!-----------------------------------------------------------------------
!     ... H2O activity
!     ... if the activity of H2O goes above 1.0, wt% can go negative
!-----------------------------------------------------------------------
               aw = ph2o_hpa / pzero_h2o
               aw = min( aw,1._r8 )
               aw = max( aw,.001_r8 )
!-----------------------------------------------------------------------
!     ... h2so4 Molality (mol/kg)
!-----------------------------------------------------------------------
               if( aw <= .05_r8 ) then
                  a1 = 12.37208932_r8
                  b1 = -0.16125516114_r8
                  c1 = -30.490657554_r8
                  d1 = -2.1133114241_r8
                  a2 = 13.455394705_r8
                  b2 = -0.1921312255_r8
                  c2 = -34.285174607_r8
                  d2 = -1.7620073078_r8
               else if( aw > .05_r8 .and. aw < .85_r8 ) then
                  a1 = 11.820654354_r8
                  b1 = -0.20786404244_r8
                  c1 = -4.807306373_r8
                  d1 = -5.1727540348_r8
                  a2 = 12.891938068_r8
                  b2 = -0.23233847708_r8
                  c2 = -6.4261237757_r8
                  d2 = -4.9005471319_r8
               else
                  a1 = -180.06541028_r8
                  b1 = -0.38601102592_r8
                  c1 = -93.317846778_r8
                  d1 = 273.88132245_r8
                  a2 = -176.95814097_r8
                  b2 = -0.36257048154_r8
                  c2 = -90.469744201_r8
                  d2 = 267.45509988_r8
               end if
!-----------------------------------------------------------------------
!     ... h2so4 mole fraction
!-----------------------------------------------------------------------
               y1       = a1*(aw**b1) + c1*aw + d1
               y2       = a2*(aw**b2) + c2*aw + d2
               m_h2so4  = y1 + ((T_limit - 190._r8)*(y2 - y1)) / 70._r8
!-----------------------------------------------------------------------
!     ... h2so4 Weight Percent
!-----------------------------------------------------------------------
               wt = 9800._r8*m_h2so4 / (98._r8*m_h2so4  + 1000._r8)
!-----------------------------------------------------------------------
!     .... Parameters for h2so4 Solution, JPL-00
!-----------------------------------------------------------------------
!     ... h2so4 Solution Density (g/cm3)
!-----------------------------------------------------------------------
	       wrk = T_limit*T_limit
               z1 =  .12364_r8  - 5.6e-7_r8*wrk
               z2 = -.02954_r8  + 1.814e-7_r8*wrk
               z3 =  2.343e-3_r8 - T_limit*1.487e-6_r8 - 1.324e-8_r8*wrk
!-----------------------------------------------------------------------
!     ... where mol_h2so4 is molality in mol/kg
!-----------------------------------------------------------------------
               den_h2so4 = 1._r8 + m_h2so4*(z1 + z2*sqrt(m_h2so4) + z3*m_h2so4)
!-----------------------------------------------------------------------
!     ... h2so4 Molarity, mol / l
!-----------------------------------------------------------------------
               molar_h2so4 = den_h2so4*wt/9.8_r8
!-----------------------------------------------------------------------
!     ... h2so4 Mole fraction
!-----------------------------------------------------------------------
               x_h2so4   = wt / (wt + ((100._r8 - wt)*98._r8/18._r8))
               term1     = .094_r8 - x_h2so4*(.61_r8 - 1.2_r8*x_h2so4)
               term2     = (8515._r8 - 10718._r8*(x_h2so4**.7_r8))*T_limiti
               H_hcl     = term1 * exp( -8.68_r8 + term2 )
               M_hcl     = H_hcl*pHCl_atm

!-----------------------------------------------------------------------
!     ... h2so4 solution viscosity
!-----------------------------------------------------------------------
               aconst    = 169.5_r8 + wt*(5.18_r8 - wt*(.0825_r8 - 3.27e-3_r8*wt))
               tzero     = 144.11_r8 + wt*(.166_r8 - wt*(.015_r8 - 2.18e-4_r8*wt))
               vis_h2so4 = aconst/(T_limit**1.43_r8) * exp( 448._r8/(T_limit - tzero) )

!-----------------------------------------------------------------------
!     ... Acid activity in molarity
!-----------------------------------------------------------------------
               term1 = 60.51_r8
               term2 = .095_r8*wt
	       wrk   = wt*wt
               term3 = .0077_r8*wrk
               term4 = 1.61e-5_r8*wt*wrk
               term5 = (1.76_r8 + 2.52e-4_r8*wrk) * T_limitsq
               term6 = -805.89_r8 + (253.05_r8*(wt**.076_r8))
               term7 = T_limitsq
               ah    = exp( term1 - term2 + term3 - term4 - term5 + term6/term7 )
	       if( ah <= 0._r8 ) then
	          write(iulog,*) 'ratecon: ah <= 0 at i,k, = ',i,k
	          write(iulog,*) 'ratecon: term1,term2,term3,term4,term5,term6,term7,wt,T_limit,ah = ', &
	                               term1,term2,term3,term4,term5,term6,term7,wt,T_limit,ah 
	       end if

	       wrk      = .25_r8*sadsulf
               rad_sulf = max( rad_sulfate(i,k),1.e-7_r8 )
!-----------------------------------------------------------------------
!     N2O5 + H2O(liq) =>  2.00*HNO3  Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
               if( n2o5vmr > small ) then
                  term0 = -25.5265_r8 - wt*(.133188_r8 - wt*(.00930846_r8 - 9.0194e-5_r8*wt))
                  term1 = 9283.76_r8 + wt*(115.345_r8 - wt*(5.19258_r8 - .0483464_r8*wt))
                  term2 = -851801._r8 - wt*(22191.2_r8 - wt*(766.916_r8 - 6.85427_r8*wt))
                  gprob_n2o5 = exp( term0 + T_limiti*(term1 + term2*T_limiti) )
                  rxt(i,k,rid_het1) = max( 0._r8,wrk*av_n2o5*gprob_n2o5 )
               end if

!-----------------------------------------------------------------------
!     ClONO2 + H2O(liq) =  HOCl + HNO3   Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
! 	... NOTE: Aerosol radius in units of cm.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     	... Radius sulfate set (from sad module)
!           Set min radius to 0.001 microns (1e-7 cm)
!           Typical radius is 0.1 microns (1e-5 cm)
!           f_cnt may go negative under if not set.
!-----------------------------------------------------------------------
                  C_cnt         = 1474._r8*T_limitsq
                  S_cnt         = .306_r8 + 24._r8*T_limiti
                  term1         = exp( -S_cnt*molar_h2so4 )
                  H_cnt         = 1.6e-6_r8 * exp( 4710._r8*T_limiti )*term1
                  D_cnt         = 5.e-8_r8*T_limit / vis_h2so4
                  k_h           = 1.22e12_r8*exp( -6200._r8*T_limiti )
                  k_h2o         = 1.95e10_r8*exp( -2800._r8*T_limiti )
                  k_hydr        = (k_h2o + k_h*ah)*aw
                  k_hcl         = 7.9e11_r8*ah*D_cnt*M_hcl
                  rdl_cnt       = sqrt( D_cnt/(k_hydr + k_hcl) )
                  term1         = 1._r8/tanh( rad_sulf/rdl_cnt )
                  term2         = rdl_cnt/rad_sulf
                  f_cnt         = term1 - term2
                  if( f_cnt > 0._r8 ) then
                     term1         = 4._r8*H_cnt*.082_r8*T_limit
                     term2         = sqrt( D_cnt*k_hydr )
                     Gamma_b_h2o   = term1*term2/C_cnt
                     term1         = sqrt( 1._r8 + k_hcl/k_hydr )
                     Gamma_cnt_rxn = f_cnt*Gamma_b_h2o*term1
                     Gamma_b_hcl   = Gamma_cnt_rxn*k_hcl/(k_hcl + k_hydr)
                     term1         = exp( -1374._r8*T_limiti )
                     Gamma_s       = 66.12_r8*H_cnt*M_hcl*term1
		     if( pHCl_atm > 0._r8 ) then
                        term1      = .612_r8*(Gamma_s*Gamma_b_hcl)* pCNT_atm/pHCl_atm
                        Fhcl       = 1._r8/(1._r8 + term1)
		     else
                        Fhcl       = 1._r8
		     end if
                     Gamma_s_prime     = Fhcl*Gamma_s
                     Gamma_b_hcl_prime = Fhcl*Gamma_b_hcl
                     term1         = Gamma_cnt_rxn*k_hydr
                     term2         = k_hcl + k_hydr
                     Gamma_b       = Gamma_b_hcl_prime + (term1/term2)
                     term1         = 1._r8 / (Gamma_s_prime + Gamma_b)
                     gprob_cnt     = 1._r8 / (1._r8 + term1)
                     term1         = Gamma_s_prime + Gamma_b_hcl_prime
                     term2         = Gamma_s_prime + Gamma_b
                     gprob_cnt_hcl = gprob_cnt * term1/term2
                     gprob_cnt_h2o = gprob_cnt - gprob_cnt_hcl
                  else
                     gprob_cnt_h2o = 0._r8
                     gprob_cnt_hcl = 0._r8
                     Fhcl          = 1._r8
                  end if
                  if( clono2vmr > small ) then
                     rxt(i,k,rid_het2) = max( 0._r8,wrk*av_clono2*gprob_cnt_h2o )
                  end if

!-----------------------------------------------------------------------
!  	... BrONO2 + H2O(liq) =  HOBr + HNO3   Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
               if( brono2vmr > small ) then
                  h1    = 29.24_r8
                  h2    = -.396_r8
                  h3    = .114_r8
                  alpha = .805_r8
                  gprob_rxn = exp( h1 + h2*wt ) + h3
                  term1     = 1._r8/alpha
                  term2     = 1._r8/gprob_rxn
                  gprob_bnt_h2o = 1._r8 / (term1 + term2)
                  rxt(i,k,rid_het3) = max( 0._r8,wrk*av_brono2*gprob_bnt_h2o )
               end if

!-----------------------------------------------------------------------
!     	... ClONO2 + HCl(liq) =  Cl2  + HNO3  Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
               if( hclvmr > small .and. clono2vmr > small ) then
                 if ( hclvmr > clono2vmr ) then
                    rxt(i,k,rid_het4) = max( 0._r8,wrk*av_clono2*gprob_cnt_hcl )*hcldeni
                 else
                    rxt(i,k,rid_het4) = max( 0._r8,wrk*av_clono2*gprob_cnt_hcl )*cntdeni
                 end if
               end if

!-----------------------------------------------------------------------
!     	... HOCl + HCl(liq) =  Cl2 + H2O   Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
               if( hclvmr > small .and. hoclvmr > small ) then
!-----------------------------------------------------------------------
!     	... Radius sulfate set (from sad module)
!           Set min radius to 0.001 microns (1e-7 cm)
!           Typical radius is 0.1 microns (1e-5 cm)
!           f_hocl may go negative under if not set.
!-----------------------------------------------------------------------
	          if( pCNT_atm > 0._r8 ) then
                     D_hocl          = 6.4e-8_r8*T_limit/vis_h2so4
                     k_hocl_hcl      = 1.25e9_r8*ah*D_hocl*M_hcl
                     C_hocl          = 2009._r8*T_limitsq
                     S_hocl          = .0776_r8 + 59.18_r8*T_limiti
                     term1           = exp( -S_hocl*molar_h2so4 )
                     H_hocl          = 1.91e-6_r8 * exp( 5862.4_r8*T_limiti )*term1
                     term1           = 4._r8*H_hocl*.082_r8*T_limit
                     term2           = sqrt( D_hocl*k_hocl_hcl )
                     Gamma_hocl_rxn  = term1*term2/C_hocl
                     rdl_hocl        = sqrt( D_hocl/k_hocl_hcl )
                     term1           = 1._r8/tanh( rad_sulf/rdl_hocl )
                     term2           = rdl_hocl/rad_sulf
                     f_hocl          = term1 - term2
                     if( f_hocl > 0._r8 ) then
                        term1           = 1._r8 / (f_hocl*Gamma_hocl_rxn*Fhcl)
                        gprob_hocl_hcl  = 1._r8 / (1._r8 + term1)
                     else
                        gprob_hocl_hcl  = 0._r8
                     end if

                     if ( hclvmr > hoclvmr ) then
                       rxt(i,k,rid_het5) = max( 0._r8,wrk*av_hocl*gprob_hocl_hcl )*hcldeni
                     else
                       rxt(i,k,rid_het5) = max( 0._r8,wrk*av_hocl*gprob_hocl_hcl )*hocldeni
                     end if

	          end if
               end if

!-----------------------------------------------------------------------
!     	... HOBr + HCl(liq) =  BrCl + H2O  Sulfate Aerosol Reaction
!-----------------------------------------------------------------------
               if( hclvmr > small .and. hobrvmr > small ) then
!-----------------------------------------------------------------------
!   	... Radius sulfate set (from sad module)
!           Set min radius to 0.001 microns (1e-7 cm)
!           Typical radius is 0.1 microns (1e-5 cm)
!           f_hobr may go negative under if not set.
!-----------------------------------------------------------------------
                  C_hobr          = 1477._r8*T_limitsq
                  D_hobr          = 9.e-9_r8
!-----------------------------------------------------------------------
!     	...  Taken from Waschewsky and Abbat
!            Dave Hanson (PC) suggested we divide this rc by eight to agee
!            with his data (Hanson, in press, 2002).
!            k1=k2*Mhcl for gamma(HOBr)
!-----------------------------------------------------------------------
                  k_wasch         = .125_r8 * exp( .542_r8*wt - 6440._r8*T_limiti + 10.3_r8)
!-----------------------------------------------------------------------
!     	... Taken from Hanson 2002.
!-----------------------------------------------------------------------
                  H_hobr          = exp( -9.86_r8 + 5427._r8*T_limiti )
                  k_dl            = 7.5e14_r8*D_hobr*2._r8                        ! or  7.5e14*D *(2nm)
!-----------------------------------------------------------------------
!  	... If k_wasch is GE than the diffusion limit...
!-----------------------------------------------------------------------
		  if( M_hcl > 0._r8 ) then
                     if( k_wasch >= k_dl ) then
                        k_hobr_hcl   = k_dl * M_hcl
                     else
                        k_hobr_hcl   = k_wasch * M_hcl
                     end if
                     term1           = 4._r8*H_hobr*.082_r8*T_limit
                     term2           = sqrt( D_hobr*k_hobr_hcl )
                     tmp             = rad_sulf/term2
                     Gamma_hobr_rxn  = term1*term2/C_hobr
                     rdl_hobr        = sqrt( D_hobr/k_hobr_hcl )
		     if( tmp < 1.e2_r8 ) then
                        term1           = 1._r8/tanh( rad_sulf/rdl_hobr )
		     else
                        term1           = 1._r8
		     end if
                     term2           = rdl_hobr/rad_sulf
                     f_hobr          = term1 - term2
                     if( f_hobr > 0._r8 ) then
                        term1            = 1._r8 / (f_hobr*Gamma_hobr_rxn)
                        gprob_hobr_hcl   = 1._r8 / (1._r8 + term1)
                     else
                         gprob_hobr_hcl  = 0._r8
                     end if

                     if ( hclvmr > hobrvmr ) then
                        rxt(i,k,rid_het6) = max( 0._r8,wrk*av_hobr*gprob_hobr_hcl )*hcldeni
                     else
                        rxt(i,k,rid_het6) = max( 0._r8,wrk*av_hobr*gprob_hobr_hcl )*hobrdeni    
                     end if           

		  end if
               end if
            end if has_sadsulf

has_sadnat : &
	    if( sadnat > 0._r8 ) then
	       wrk = .25_r8*sadnat
!-----------------------------------------------------------------------
!     	... N2O5 + H2O(s) => 2HNO3  NAT Aerosol Reaction
!-----------------------------------------------------------------------
               if( n2o5vmr > small ) then
!-----------------------------------------------------------------------
!     ... gprob based on JPL06-2 for NAT.
!         also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
!                 gprob_tot     = 4.e-4
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het7)  = wrk*av_n2o5*4.e-4_r8
	       end if

!-----------------------------------------------------------------------
!     ClONO2 + H2O(s) => HNO3 + HOCl  NAT Aerosol Reaction
!-----------------------------------------------------------------------
               if( clono2vmr > small ) then
!-----------------------------------------------------------------------
!     ... gprob based on JPL06-2 for NAT.
!         also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
!                 gprob_tot    = 0.004
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het8) = wrk*av_clono2*4.0e-3_r8
   	       end if

!-----------------------------------------------------------------------
!     	... ClONO2 + HCl(s) => HNO3 + Cl2, NAT Aerosol Reaction
!-----------------------------------------------------------------------
               if( hclvmr > small ) then
                  if( clono2vmr > small ) then
!-----------------------------------------------------------------------
!     ... gprob based on JPL06-2 for NAT.
!         also see Hanson and Ravi, JPC, 96, 2682-2691, 1992.
!                 gprob_tot   = 0.2
!-----------------------------------------------------------------------
                     if ( hclvmr > clono2vmr ) then
                        rxt(i,k,rid_het9) = wrk*av_clono2*0.2_r8*hcldeni
                     else
                        rxt(i,k,rid_het9) = wrk*av_clono2*0.2_r8*cntdeni  
                     end if
                  end if

!-----------------------------------------------------------------------
!     	... HOCl + HCl(s) => H2O + Cl2  NAT Aerosol Reaction
!-----------------------------------------------------------------------
                  if( hoclvmr > small ) then
!-----------------------------------------------------------------------
!     ... gprob based on JPL06-2 for NAT.
!         see Hanson and Ravi, JPC, 96, 2682-2691, 1992.
!         and      Abbatt and Molina, GRL, 19, 461-464, 1992.
!                 gprob_tot   = 0.1
!-----------------------------------------------------------------------
                     if ( hclvmr > hoclvmr ) then
                        rxt(i,k,rid_het10) = wrk*av_hocl*0.1_r8*hcldeni
                     else
                        rxt(i,k,rid_het10) = wrk*av_hocl*0.1_r8*hocldeni
                     end if
                  end if
               end if

!-----------------------------------------------------------------------
!     	... BrONO2 + H2O(s) => HOBr + HNO3  NAT Aerosol Reaction
!-----------------------------------------------------------------------
               if( brono2vmr > small ) then
!-----------------------------------------------------------------------
!       ... Personel Communication, 11/4/99, David Hanson
!                 gprob_tot   = 0.3
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het11) = wrk*av_brono2*0.3_r8
               end if
            end if has_sadnat

has_sadice : &
	    if( sadice > 0._r8 ) then
	       wrk = .25_r8*sadice
!-----------------------------------------------------------------------
!     N2O5 + H2O(s) => 2HNO3  ICE Aerosol Reaction
!-----------------------------------------------------------------------
               if( n2o5vmr > small ) then
!-----------------------------------------------------------------------
!       ... gprob based on JPL06-2 for ICE.
!           also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
!                 gprob_tot    = .02
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het12) = wrk*av_n2o5*0.02_r8
 	       end if
!-----------------------------------------------------------------------
!     	... ClONO2 + H2O(s) => HNO3 + HOCl  ICE Aerosol Reaction
!-----------------------------------------------------------------------
               if( clono2vmr > small ) then
!-----------------------------------------------------------------------
!     	... gprob based on JPL06-2 for ICE.
!     	    also see Hanson and Ravi, JGR, 96, 17307-17314, 1991.
!                 gprob_tot    = .3
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het13) = wrk*av_clono2*0.3_r8
	       end if

!-----------------------------------------------------------------------
!     	... BrONO2 + H2O(s) => HNO3 + HOBr  ICE Aerosol Reaction
!-----------------------------------------------------------------------
               if( brono2vmr > small ) then
!-----------------------------------------------------------------------
!     	... gprob based on JPL06-2 for ICE.
!           also see Hanson and Ravi, JPC, 97, 2802-2803, 1993.
!           could be as high as 1.0
!                 gprob_tot    = .3
!-----------------------------------------------------------------------
                  rxt(i,k,rid_het14) = wrk*av_brono2*0.3_r8
      	       end if

!-----------------------------------------------------------------------
!     ClONO2 + HCl(s) => HNO3 + Cl2, ICE Aerosol Reaction
!-----------------------------------------------------------------------
               if( hclvmr > small ) then
                  if( clono2vmr > small ) then
!-----------------------------------------------------------------------
!       ... gprob based on JPL06-2 for ICE.
!           also see Hanson and Ravi, GRL, 15, 17-20, 1988.
!           also see Lue et al.,
!                 gprob_tot    = .3
!-----------------------------------------------------------------------
                     if ( hclvmr > clono2vmr ) then
                        rxt(i,k,rid_het15) = wrk*av_clono2*0.3_r8*hcldeni
                     else
                        rxt(i,k,rid_het15) = wrk*av_clono2*0.3_r8*cntdeni
                     end if

                  end if
!
!-----------------------------------------------------------------------
!     	... HOCl + HCl(s) => H2O + Cl2, ICE Aerosol Reaction
!-----------------------------------------------------------------------
                  if( hoclvmr > small .and. hclvmr > small ) then
!-----------------------------------------------------------------------
!       ... gprob based on JPL06-2 for ICE.
!           also see Hanson and Ravi, JPC, 96, 2682-2691, 1992.
!           also see Abbatt and Molina, GRL, 19, 461-464, 1992.
!                 gprob_tot   = .2
!-----------------------------------------------------------------------
                     if ( hclvmr > hoclvmr ) then
                        rxt(i,k,rid_het16) = wrk*av_hocl*0.2_r8*hcldeni
                     else
                        rxt(i,k,rid_het16) = wrk*av_hocl*0.2_r8*hocldeni
                     end if

                  end if

!-----------------------------------------------------------------------
!     HOBr + HCl(s) => H2O + BrCl, ICE Aerosol Reaction
!-----------------------------------------------------------------------
                  if( hobrvmr > small .and. hclvmr > small ) then
!-----------------------------------------------------------------------
!       ... gprob based on JPL06-2 for ICE.
!           Abbatt GRL, 21, 665-668, 1994.
!                    gprob_tot   = .3
!-----------------------------------------------------------------------
                    if ( hclvmr > hobrvmr ) then
                       rxt(i,k,rid_het17) = wrk*av_hobr*0.3_r8*hcldeni
                    else
                       rxt(i,k,rid_het17) = wrk*av_hobr*0.3_r8*hobrdeni
                    end if
                  end if
	       end if
            end if has_sadice
         end do column_loop
      end do Level_loop

      end subroutine ratecon_sfstrat

      end module mo_strato_rates
