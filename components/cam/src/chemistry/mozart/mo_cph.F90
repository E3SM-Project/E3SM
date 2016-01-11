
      module mo_cph
!---------------------------------------------------------------------
!	... compute chemical potential heating
!---------------------------------------------------------------------

      use shr_kind_mod,  only : r8 => shr_kind_r8

      implicit none

      save

      integer , parameter :: ncph    = 41
      real(r8), parameter :: secpday = 86400._r8
      real(r8), parameter :: daypsec = 1._r8/secpday
!==============================================================
!... Doug Kinnison, dkin@ucar.edu
!
!... Enthalpy Data are taken from Atkinson et al., 
!    Evaluated kinetic and photochemical data for atmospheric
!    chemistry: Volume I, Atmos. Chem. Phys., 4, 1461-1738.
!... Heats of formation at 0K.
!... Units: kJ mol-1
!
!... Exception to the Atkinson et al. reference  (@0K unless noted)
!    (4), (5), (8), (9), (10, (11), (14), (15), (27), (28) at 298K  
!    (7)  h + o2 -> oh + o2 is multiplied by 0.6 (Mlynczak) to represent
!         AG loss of excited OH.
!    (25) n2d + o2 -> no + o1d taken from Roble, UMLT, Johnson and Killeen.
!    (26) n2d + o  -> n  + o   taken from Roble, UMLT, Johnson and Killeen.
!    (30-41) Taken from Roble, UMLT, Johnson and Killeen Ed., Geophys. Mono. 87
!==============================================================
      real(r8), parameter :: exotherm(ncph) = (/ &
         392.19_r8, &                       !(1)  o + o3 -> 2o2
         493.58_r8, &                       !(2)  o + o + m -> o2 + m
         67.67_r8,  &                       !(3)  o + oh -> h + o2
         226.58_r8, &                       !(4)  o + ho2 -> oh + o2
         203.40_r8, &                       !(5)  h + o2 + m -> ho2 + m
         101.39_r8, &                       !(6)  o + o2 + m -> o3 + m
         194.71_r8, &                       !(7)  h + o3 -> oh + o2 (*0.6)
         34.47_r8,  &                       !(8)  ho2 + no -> no2 + oh
         120.10_r8, &                       !(9)  ho2 + o3 -> oh + 2o2
         165.51_r8, &                       !(10) ho2 + ho2 -> h2o2 + o2
         165.30_r8, &                       !(11) oh + o3 -> ho2 + o2
         199.17_r8, &                       !(12) no + o3 -> no2 + o2
         193.02_r8, &                       !(13) no2 + o -> no + o2
         293.62_r8, &                       !(14) oh + ho2 -> h2o + o2
         232.59_r8, &                       !(15) h + ho2 -> h2 + o2
         32.91_r8,  &                       !(16) o1d + o2 -> o + o2_1s
         189.91_r8, &                       !(17) o1d + n2 -> o + n2
         62.60_r8,  &                       !(18) o2_1s + o -> o2_1d + o
         62.60_r8,  &                       !(19) o2_1s + o2 -> o2_1d + o2
         62.60_r8,  &                       !(20) o2_1s + n2 -> o2_1d + n2
         62.60_r8,  &                       !(21) o2_1s + o3 -> o2_1d + o3
         94.30_r8,  &                       !(22) o2_1d + o -> o2 + o
         94.30_r8,  &                       !(23) o2_1d + o2 -> 2o2
         94.30_r8,  &                       !(24) o2_1d + n2 -> o2 + n2
         177.51_r8, &                       !(25) n2d + o2 -> no + o1d
         229.61_r8, &                       !(26) n2d + o -> n + o
         133.75_r8, &                       !(27) n + o2 -> no + o
         313.75_r8, &                       !(28) n + no -> n2 + o
         189.81_r8, &                       !(29) o1d + o2 -> o + o2
         150.11_r8, &                       !(30) Op + o2 -> O2p + o
         105.04_r8, &                       !(31) Op + n2 -> NOp + n
         67.53_r8,  &                       !(32) N2p + o -> NOp + n2d
         406.16_r8, &                       !(33) O2p + n -> NOp + o
         271.38_r8, &                       !(34) O2p + no -> NOp + o2
         239.84_r8, &                       !(35) Np + o2 -> O2p + n
         646.28_r8, &                       !(36) Np + o2 -> NOp + o
         95.55_r8,  &                       !(37) Np + o -> Op + n
         339.59_r8, &                       !(38) N2p + o2 -> O2p + n2
         82.389_r8, &                       !(39) NOp + e -> 
         508.95_r8, &                       !(40) O2p + e -> 
         354.83_r8  &                       !(41) N2p + e -> 
         /)

      private
      public :: cph, init_cph

      integer :: rid_cph1,rid_cph2,rid_cph3,rid_cph4,rid_cph5,rid_cph6,rid_cph7,rid_cph8, &
                 rid_cph9,rid_cph10,rid_cph11,rid_cph12,rid_cph13,rid_cph14,rid_cph15,&
                 rid_cph16,rid_cph17,rid_cph18,rid_cph19,rid_cph20,rid_cph21,rid_cph22,&
                 rid_cph23,rid_cph24,rid_cph25,rid_cph26,rid_cph27,rid_cph28,rid_cph29
      integer :: rid_ion1,rid_ion2,rid_ion3,rid_ion4,rid_ion5,rid_ion6,rid_ion7,rid_ion8,rid_ion9
      integer :: rid_elec1,rid_elec2,rid_elec3

      integer :: id_o3, id_o, id_h, id_ho2, id_oh, id_no, id_no2, id_o1d, id_o2_1s, id_o2_1d, id_n2d, id_n, &
                 id_op, id_np, id_nop, id_n2p, id_o2p, id_o2, id_e

      logical :: has_cph

      contains

       subroutine init_cph

          use mo_chem_utls, only : get_rxt_ndx, get_spc_ndx
          use cam_history,  only : addfld
          use ppgrid,       only : pver

          implicit none

          integer :: ids(60) = -1

          rid_cph1  = get_rxt_ndx( 'cph_O_O3' )
          rid_cph2  = get_rxt_ndx( 'usr_O_O' )
          rid_cph3  = get_rxt_ndx( 'cph_OH_O' )
          rid_cph4  = get_rxt_ndx( 'cph_HO2_O' )
          rid_cph5  = get_rxt_ndx( 'cph_H_O2' )
          rid_cph6  = get_rxt_ndx( 'usr_O_O2' )
          rid_cph7  = get_rxt_ndx( 'cph_H_O3' )
          rid_cph8  = get_rxt_ndx( 'cph_NO_HO2' )
          rid_cph9  = get_rxt_ndx( 'cph_HO2_O3' )
          rid_cph10 = get_rxt_ndx( 'usr_HO2_HO2' )
          rid_cph11 = get_rxt_ndx( 'cph_OH_O3' )
          rid_cph12 = get_rxt_ndx( 'cph_NO_O3' )
          rid_cph13 = get_rxt_ndx( 'cph_NO2_O' )
          rid_cph14 = get_rxt_ndx( 'cph_OH_HO2' )
          rid_cph15 = get_rxt_ndx( 'cph_H_HO2' )
          rid_cph16 = get_rxt_ndx( 'cph_O1D_O2' )
          rid_cph17 = get_rxt_ndx( 'cph_O1D_N2' )
          rid_cph18 = get_rxt_ndx( 'cph_O2_1S_O' )
          rid_cph19 = get_rxt_ndx( 'cph_O2_1S_O2' )
          rid_cph20 = get_rxt_ndx( 'cph_O2_1S_N2' )
          rid_cph21 = get_rxt_ndx( 'cph_O2_1S_O3' )
          rid_cph22 = get_rxt_ndx( 'cph_O2_1D_O' )
          rid_cph23 = get_rxt_ndx( 'cph_O2_1D_O2' )
          rid_cph24 = get_rxt_ndx( 'cph_O2_1D_N2' )
          rid_cph25 = get_rxt_ndx( 'cph_N2D_O2' )
          rid_cph26 = get_rxt_ndx( 'cph_N2D_O' )
          rid_cph27 = get_rxt_ndx( 'cph_N_O2' )
          rid_cph28 = get_rxt_ndx( 'cph_N_NO' )
          rid_cph29 = get_rxt_ndx( 'cph_O1D_O2b' )
          rid_ion1  = get_rxt_ndx( 'ion_Op_O2' )
          rid_ion2  = get_rxt_ndx( 'ion_Op_N2' )
          rid_ion3  = get_rxt_ndx( 'ion_N2p_Oa' )
          rid_ion4  = get_rxt_ndx( 'ion_O2p_N' )
          rid_ion5  = get_rxt_ndx( 'ion_O2p_NO' )
          rid_ion6  = get_rxt_ndx( 'ion_Np_O2a' )
          rid_ion7  = get_rxt_ndx( 'ion_Np_O2b' )
          rid_ion8  = get_rxt_ndx( 'ion_Np_O' )
          rid_ion9  = get_rxt_ndx( 'ion_N2p_O2' )
          rid_elec1 = get_rxt_ndx( 'elec1' )
          rid_elec2 = get_rxt_ndx( 'elec2' )
          rid_elec3 = get_rxt_ndx( 'elec3' )
          id_o3     = get_spc_ndx( 'O3' )
          id_o      = get_spc_ndx( 'O' )
          id_h      = get_spc_ndx( 'H' )
          id_ho2    = get_spc_ndx( 'HO2' )
          id_oh     = get_spc_ndx( 'OH' )
          id_no     = get_spc_ndx( 'NO' )
          id_no2    = get_spc_ndx( 'NO2' )
          id_o1d    = get_spc_ndx( 'O1D' )
          id_o2     = get_spc_ndx( 'O2' )
          id_o2_1s  = get_spc_ndx( 'O2_1S' )
          id_o2_1d  = get_spc_ndx( 'O2_1D' )
          id_n2d    = get_spc_ndx( 'N2D' )
          id_n      = get_spc_ndx( 'N' )
          id_op     = get_spc_ndx( 'Op' )
          id_nop    = get_spc_ndx( 'NOp' )
          id_n2p    = get_spc_ndx( 'N2p' )
          id_np     = get_spc_ndx( 'Np' )
          id_o2p    = get_spc_ndx( 'O2p' )
          id_e      = get_spc_ndx( 'e' )

          ids(:) = (/ rid_cph1,rid_cph2,rid_cph3,rid_cph4,rid_cph5,rid_cph6,rid_cph7,rid_cph8, &
                 rid_cph9,rid_cph10,rid_cph11,rid_cph12,rid_cph13,rid_cph14,rid_cph15,&
                 rid_cph16,rid_cph17,rid_cph18,rid_cph19,rid_cph20,rid_cph21,rid_cph22,&
                 rid_cph23,rid_cph24,rid_cph25,rid_cph26,rid_cph27,rid_cph28,rid_cph29,&
                 rid_ion1,rid_ion2,rid_ion3,rid_ion4,rid_ion5,rid_ion6,rid_ion7,rid_ion8,rid_ion9, &
                 rid_elec1,rid_elec2,rid_elec3,&
                 id_o3, id_o, id_h, id_ho2, id_oh, id_no, id_nop, id_no2, id_o1d, id_o2_1s, id_o2_1d, id_n2d, id_n, &
                 id_op, id_n2p, id_np, id_o2p, id_o2, id_e  /)

          has_cph = any( ids(:) > 0 )

          if (.not.has_cph) return

          call addfld( 'CPH1', (/ 'lev' /), 'I',  'K/s', 'O + O3 -> 2*O2 + 93.46 kcal/mol chem pot heating rate' )
          call addfld( 'CPH2', (/ 'lev' /), 'I',  'K/s', 'O + O + M -> O2 + M + 119.12 kcal/mol chem pot heating rate' )
          call addfld( 'CPH3', (/ 'lev' /), 'I',  'K/s', 'O + OH -> H + O2 + 16.35 kcal/mol chem pot heating rate' )
          call addfld( 'CPH4', (/ 'lev' /), 'I',  'K/s', 'O + HO2 -> OH + O2 + 53.97 kcal/mol chem pot heating rate' )
          call addfld( 'CPH5', (/ 'lev' /), 'I',  'K/s', 'H + O2 + M -> HO2 + M + 48.8 kcal/mol chem pot heating rate' )
          call addfld( 'CPH6', (/ 'lev' /), 'I',  'K/s', 'O + O2 + M -> O3 + M + 25.66 kcal/mol chem pot heating rate' )
          call addfld( 'CPH7', (/ 'lev' /), 'I',  'K/s', 'H + O3 -> OH + O2 + 46.27 kcal/mol chem pot heating rate' )
          call addfld( 'CPH8', (/ 'lev' /), 'I',  'K/s', 'HO2 + NO -> NO2 + OH + 8.06 kcal/mol chem pot heating rate' )
          call addfld( 'CPH9', (/ 'lev' /), 'I',  'K/s', 'HO2 + O3 -> 2*O2 + OH + 28.31 kcal/mol chem pot heating rate' )
          call addfld( 'CPH10', (/ 'lev' /), 'I', 'K/s', 'HO2 + HO2 -> H2O2 + O2 + 39.08 kcal/mol chem pot heating rate' )
          call addfld( 'CPH11', (/ 'lev' /), 'I', 'K/s', 'OH + O3 -> HO2 + O2 + 39.49 kcal/mol chem pot heating rate' )
          call addfld( 'CPH12', (/ 'lev' /), 'I', 'K/s', 'NO + O3 -> NO2 + O2 + 47.55 kcal/mol chem pot heating rate' )
          call addfld( 'CPH13', (/ 'lev' /), 'I', 'K/s', 'NO2 + O -> NO + O2 + 45.91 kcal/mol chem pot heating rate' )
          call addfld( 'CPH14', (/ 'lev' /), 'I', 'K/s', 'OH + HO2 -> H2O + O2 + 69.99 kcal/mol chem pot heating rate' )
          call addfld( 'CPH15', (/ 'lev' /), 'I', 'K/s', 'H + HO2 -> H2 + O2 + 55.4 kcal/mol chem pot heating rate' )
          call addfld( 'CPH16', (/ 'lev' /), 'I', 'K/s', 'O1D + O2 -> O + O2_1S + 7.63 kcal/mol chem pot heating rate' )
          call addfld( 'CPH17', (/ 'lev' /), 'I', 'K/s', 'O1D + N2 -> O + N2 + 45.14 kcal/mol chem pot heating rate' )
          call addfld( 'CPH18', (/ 'lev' /), 'I', 'K/s', 'O2_1S + O -> O2_1D + O + 14.97 kcal/mol chem pot heating rate' )
          call addfld( 'CPH19', (/ 'lev' /), 'I', 'K/s', 'O2_1S + O2 -> O2_1D + O2 + 14.97 kcal/mol chem pot heating rate' )
          call addfld( 'CPH20', (/ 'lev' /), 'I', 'K/s', 'O2_1S + N2 -> O2_1D + N2 + 14.97 kcal/mol chem pot heating rate' )
          call addfld( 'CPH21', (/ 'lev' /), 'I', 'K/s', 'O2_1S + O3 -> O2_1D + O3 + 14.97 kcal/mol chem pot heating rate' )
          call addfld( 'CPH22', (/ 'lev' /), 'I', 'K/s', 'O2_1D + O -> O2 + O + 22.54 kcal/mol chem pot heating rate' )
          call addfld( 'CPH23', (/ 'lev' /), 'I', 'K/s', 'O2_1D + O2 -> 2*O2 + 22.54 kcal/mol chem pot heating rate' )
          call addfld( 'CPH24', (/ 'lev' /), 'I', 'K/s', 'O2_1D + N2 -> O2 + N2 + 22.54 kcal/mol chem pot heating rate' )
          call addfld( 'CPH25', (/ 'lev' /), 'I', 'K/s', 'N2D + O2 -> NO + O1D + 42.43 kcal/mol chem pot heating rate' )
          call addfld( 'CPH26', (/ 'lev' /), 'I', 'K/s', 'N2D + O -> N + O + 54.88 kcal/mol chem pot heating rate' )
          call addfld( 'CPH27', (/ 'lev' /), 'I', 'K/s', 'N + O2 -> NO + O + 32.28 kcal/mol chem pot heating rate' )
          call addfld( 'CPH28', (/ 'lev' /), 'I', 'K/s', 'N + NO -> N2 + O + 61.8 kcal/mol chem pot heating rate' )
          call addfld( 'CPH29', (/ 'lev' /), 'I', 'K/s', 'O1D + O2 -> O + O2 + 45.14 kcal/mol chem pot heating rate' )
          call addfld( 'CPH30', (/ 'lev' /), 'I', 'K/s', 'Op + O2 -> O2p + O + 35.878 kcal/mol chem pot heating rate' )
          call addfld( 'CPH31', (/ 'lev' /), 'I', 'K/s', 'Op + N2 -> NOp + N + 25.105 kcal/mol chem pot heating rate' )
          call addfld( 'CPH32', (/ 'lev' /), 'I', 'K/s', 'N2p + O -> NOp + N2D + 16.140 kcal/mol chem pot heating rate' )
          call addfld( 'CPH33', (/ 'lev' /), 'I', 'K/s', 'O2p + N -> NOp + O + 97.073 kcal/mol chem pot heating rate' )
          call addfld( 'CPH34', (/ 'lev' /), 'I', 'K/s', 'O2p + NO -> NOp + O2 + 64.862 kcal/mol chem pot heating rate' )
          call addfld( 'CPH35', (/ 'lev' /), 'I', 'K/s', 'Np + O2 -> O2p + N + 57.322 kcal/mol chem pot heating rate' )
          call addfld( 'CPH36', (/ 'lev' /), 'I', 'K/s', 'Np + O2 -> NOp + O + 154.46 kcal/mol chem pot heating rate' )
          call addfld( 'CPH37', (/ 'lev' /), 'I', 'K/s', 'Np + O -> Op + N + 22.597 kcal/mol chem pot heating rate' )
          call addfld( 'CPH38', (/ 'lev' /), 'I', 'K/s', 'N2p + O2 -> O2p + N2 + 81.164 kcal/mol chem pot heating rate' )
          call addfld( 'CPH39', (/ 'lev' /), 'I', 'K/s', 'NOp + e -> 19.691 kcal/mol chem pot heating rate' )
          call addfld( 'CPH40', (/ 'lev' /), 'I', 'K/s', 'O2p + e -> 121.65 kcal/mol chem pot heating rate' )
          call addfld( 'CPH41', (/ 'lev' /), 'I', 'K/s', 'N2p + e -> 84.807 kcal/mol chem pot heating rate' )
          call addfld( 'QCP', (/ 'lev' /), 'I',   'K/s', 'chem pot heating rate' )

       end subroutine init_cph

      subroutine cph( cph_tot, vmr, rxt, cp, mbar, &
                      kbot, ncol, lchnk )
!-----------------------------------------------------------------------
!      	... forms the chemical potential heating rates
!-----------------------------------------------------------------------

      use chem_mods,     only : gas_pcnst, rxntot
      use ppgrid,        only : pver
      use cam_history,   only : outfld

      implicit none

!-----------------------------------------------------------------------
!     	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)   ::  ncol                                ! columns in chunck
      integer, intent(in)   ::  lchnk                               ! chunk index
      integer, intent(in)   ::  kbot                                ! bottom vert index
      real(r8), intent(in)  ::  rxt(ncol,pver,rxntot)               ! rxt rates (1/cm^3/s)
      real(r8), intent(in)  ::  cp(ncol,pver)                       ! specific heat capacity (J/K/kg)
      real(r8), intent(in)  ::  mbar(ncol,pver)                     ! atm mean mass (g/mole)
      real(r8), intent(in)  ::  vmr(ncol,pver,gas_pcnst)            ! concentrations (mol/mol)
      real(r8), intent(out) ::  cph_tot(ncol,pver)                  ! total heating (K/s)

!-----------------------------------------------------------------------
!     	... local variables
!-----------------------------------------------------------------------
      integer  ::  i, k
      real(r8) ::  tmp(ncol)
      real(r8) ::  cph_rate(ncol,pver,ncph)

      if (.not.has_cph) return
       
      cph_rate(:,:,:) = 0._r8
      do k = 1,kbot
         tmp(:)           = 1._r8 / (1.e-6_r8*cp(:,k)*mbar(:,k))
         if (rid_cph1 >0) cph_rate(:,k,1)  = tmp(:)*rxt(:,k,rid_cph1 )*vmr(:,k,id_o)*vmr(:,k,id_o3)*exotherm(1)
         if (rid_cph2 >0) cph_rate(:,k,2)  = tmp(:)*rxt(:,k,rid_cph2 )*vmr(:,k,id_o)*vmr(:,k,id_o)*exotherm(2)
         if (rid_cph3 >0) cph_rate(:,k,3)  = tmp(:)*rxt(:,k,rid_cph3 )*vmr(:,k,id_o)*vmr(:,k,id_oh)*exotherm(3)
         if (rid_cph4 >0) cph_rate(:,k,4)  = tmp(:)*rxt(:,k,rid_cph4 )*vmr(:,k,id_o)*vmr(:,k,id_ho2)*exotherm(4)
         if (rid_cph5 >0) cph_rate(:,k,5)  = tmp(:)*rxt(:,k,rid_cph5 )*vmr(:,k,id_h)*vmr(:,k,id_o2)*exotherm(5)
         if (rid_cph6 >0) cph_rate(:,k,6)  = tmp(:)*rxt(:,k,rid_cph6 )*vmr(:,k,id_o)*vmr(:,k,id_o2)*exotherm(6)
         if (rid_cph7 >0) cph_rate(:,k,7)  = tmp(:)*rxt(:,k,rid_cph7 )*vmr(:,k,id_h)*vmr(:,k,id_o3)*exotherm(7)
         if (rid_cph8 >0) cph_rate(:,k,8)  = tmp(:)*rxt(:,k,rid_cph8 )*vmr(:,k,id_ho2)*vmr(:,k,id_no)*exotherm(8)
         if (rid_cph9 >0) cph_rate(:,k,9)  = tmp(:)*rxt(:,k,rid_cph9 )*vmr(:,k,id_ho2)*vmr(:,k,id_o3)*exotherm(9)
         if (rid_cph10>0) cph_rate(:,k,10) = tmp(:)*rxt(:,k,rid_cph10)*vmr(:,k,id_ho2)*vmr(:,k,id_ho2)*exotherm(10)
         if (rid_cph11>0) cph_rate(:,k,11) = tmp(:)*rxt(:,k,rid_cph11)*vmr(:,k,id_oh)*vmr(:,k,id_o3)*exotherm(11)
         if (rid_cph12>0) cph_rate(:,k,12) = tmp(:)*rxt(:,k,rid_cph12)*vmr(:,k,id_no)*vmr(:,k,id_o3)*exotherm(12)
         if (rid_cph13>0) cph_rate(:,k,13) = tmp(:)*rxt(:,k,rid_cph13)*vmr(:,k,id_no2)*vmr(:,k,id_o)*exotherm(13)
         if (rid_cph14>0) cph_rate(:,k,14) = tmp(:)*rxt(:,k,rid_cph14)*vmr(:,k,id_oh)*vmr(:,k,id_ho2)*exotherm(14)
         if (rid_cph15>0) cph_rate(:,k,15) = tmp(:)*rxt(:,k,rid_cph15)*vmr(:,k,id_h)*vmr(:,k,id_ho2)*exotherm(15)
         if (rid_cph16>0) cph_rate(:,k,16) = tmp(:)*rxt(:,k,rid_cph16)*vmr(:,k,id_o1d)*vmr(:,k,id_o2)*exotherm(16)
         if (rid_cph17>0) cph_rate(:,k,17) = tmp(:)*rxt(:,k,rid_cph17)*vmr(:,k,id_o1d)*exotherm(17)
         if (rid_cph18>0) cph_rate(:,k,18) = tmp(:)*rxt(:,k,rid_cph18)*vmr(:,k,id_o2_1s)*vmr(:,k,id_o)*exotherm(18)
         if (rid_cph19>0) cph_rate(:,k,19) = tmp(:)*rxt(:,k,rid_cph19)*vmr(:,k,id_o2_1s)*vmr(:,k,id_o2)*exotherm(19)
         if (rid_cph20>0) cph_rate(:,k,20) = tmp(:)*rxt(:,k,rid_cph20)*vmr(:,k,id_o2_1s)*exotherm(20)
         if (rid_cph21>0) cph_rate(:,k,21) = tmp(:)*rxt(:,k,rid_cph21)*vmr(:,k,id_o2_1s)*vmr(:,k,id_o3)*exotherm(21)
         if (rid_cph22>0) cph_rate(:,k,22) = tmp(:)*rxt(:,k,rid_cph22)*vmr(:,k,id_o2_1d)*vmr(:,k,id_o)*exotherm(22)
         if (rid_cph23>0) cph_rate(:,k,23) = tmp(:)*rxt(:,k,rid_cph23)*vmr(:,k,id_o2_1d)*vmr(:,k,id_o2)*exotherm(23)
         if (rid_cph24>0) cph_rate(:,k,24) = tmp(:)*rxt(:,k,rid_cph24)*vmr(:,k,id_o2_1d)*exotherm(24)
         if (rid_cph25>0) cph_rate(:,k,25) = tmp(:)*rxt(:,k,rid_cph25)*vmr(:,k,id_n2d)*vmr(:,k,id_o2)*exotherm(25)
         if (rid_cph26>0) cph_rate(:,k,26) = tmp(:)*rxt(:,k,rid_cph26)*vmr(:,k,id_n2d)*vmr(:,k,id_o)*exotherm(26)
         if (rid_cph27>0) cph_rate(:,k,27) = tmp(:)*rxt(:,k,rid_cph27)*vmr(:,k,id_n)*vmr(:,k,id_o2)*exotherm(27)
         if (rid_cph28>0) cph_rate(:,k,28) = tmp(:)*rxt(:,k,rid_cph28)*vmr(:,k,id_n)*vmr(:,k,id_no)*exotherm(28)
         if (rid_cph29>0) cph_rate(:,k,29) = tmp(:)*rxt(:,k,rid_cph29)*vmr(:,k,id_o1d)*vmr(:,k,id_o2)*exotherm(29)
         if (rid_ion1 >0) cph_rate(:,k,30) = tmp(:)*rxt(:,k,rid_ion1 )*vmr(:,k,id_op)*vmr(:,k,id_o2)*exotherm(30)
         if (rid_ion2 >0) cph_rate(:,k,31) = tmp(:)*rxt(:,k,rid_ion2 )*vmr(:,k,id_op)*exotherm(31)
         if (rid_ion3 >0) cph_rate(:,k,32) = tmp(:)*rxt(:,k,rid_ion3 )*vmr(:,k,id_n2p)*vmr(:,k,id_o)*exotherm(32)
         if (rid_ion4 >0) cph_rate(:,k,33) = tmp(:)*rxt(:,k,rid_ion4 )*vmr(:,k,id_o2p)*vmr(:,k,id_n)*exotherm(33)
         if (rid_ion5 >0) cph_rate(:,k,34) = tmp(:)*rxt(:,k,rid_ion5 )*vmr(:,k,id_o2p)*vmr(:,k,id_no)*exotherm(34)
         if (rid_ion6 >0) cph_rate(:,k,35) = tmp(:)*rxt(:,k,rid_ion6 )*vmr(:,k,id_o2)*vmr(:,k,id_np)*exotherm(35)
         if (rid_ion7 >0) cph_rate(:,k,36) = tmp(:)*rxt(:,k,rid_ion7 )*vmr(:,k,id_o2)*vmr(:,k,id_np)*exotherm(36)
         if (rid_ion8 >0) cph_rate(:,k,37) = tmp(:)*rxt(:,k,rid_ion8 )*vmr(:,k,id_np)*vmr(:,k,id_o)*exotherm(37)
         if (rid_ion9 >0) cph_rate(:,k,38) = tmp(:)*rxt(:,k,rid_ion9 )*vmr(:,k,id_n2p)*vmr(:,k,id_o2)*exotherm(38)
         if (rid_elec1>0) cph_rate(:,k,39) = tmp(:)*rxt(:,k,rid_elec1)*vmr(:,k,id_e)*vmr(:,k,id_nop)*exotherm(39)
         if (rid_elec2>0) cph_rate(:,k,40) = tmp(:)*rxt(:,k,rid_elec2)*vmr(:,k,id_e)*vmr(:,k,id_o2p)*exotherm(40)
         if (rid_elec3>0) cph_rate(:,k,41) = tmp(:)*rxt(:,k,rid_elec3)*vmr(:,k,id_e)*vmr(:,k,id_n2p)*exotherm(41)
      end do

      cph_tot(:,:) = 0._r8
      do k = 1,kbot
         do i = 1,ncol
            cph_tot(i,k) = sum( cph_rate(i,k,:) )
         end do
      end do

!-----------------------------------------------------------------------
!     	... output the rates
!-----------------------------------------------------------------------
      call outfld( 'CPH1', cph_rate(:,:,1), ncol, lchnk )
      call outfld( 'CPH2', cph_rate(:,:,2), ncol, lchnk )
      call outfld( 'CPH3', cph_rate(:,:,3), ncol, lchnk )
      call outfld( 'CPH4', cph_rate(:,:,4), ncol, lchnk )
      call outfld( 'CPH5', cph_rate(:,:,5), ncol, lchnk )
      call outfld( 'CPH6', cph_rate(:,:,6), ncol, lchnk )
      call outfld( 'CPH7', cph_rate(:,:,7), ncol, lchnk )
      call outfld( 'CPH8', cph_rate(:,:,8), ncol, lchnk )
      call outfld( 'CPH9', cph_rate(:,:,9), ncol, lchnk )
      call outfld( 'CPH10', cph_rate(:,:,10), ncol, lchnk )
      call outfld( 'CPH11', cph_rate(:,:,11), ncol, lchnk )
      call outfld( 'CPH12', cph_rate(:,:,12), ncol, lchnk )
      call outfld( 'CPH13', cph_rate(:,:,13), ncol, lchnk )
      call outfld( 'CPH14', cph_rate(:,:,14), ncol, lchnk )
      call outfld( 'CPH15', cph_rate(:,:,15), ncol, lchnk )
      call outfld( 'CPH16', cph_rate(:,:,16), ncol, lchnk )
      call outfld( 'CPH17', cph_rate(:,:,17), ncol, lchnk )
      call outfld( 'CPH18', cph_rate(:,:,18), ncol, lchnk )
      call outfld( 'CPH19', cph_rate(:,:,19), ncol, lchnk )
      call outfld( 'CPH20', cph_rate(:,:,20), ncol, lchnk )
      call outfld( 'CPH21', cph_rate(:,:,21), ncol, lchnk )
      call outfld( 'CPH22', cph_rate(:,:,22), ncol, lchnk )
      call outfld( 'CPH23', cph_rate(:,:,23), ncol, lchnk )
      call outfld( 'CPH24', cph_rate(:,:,24), ncol, lchnk )
      call outfld( 'CPH25', cph_rate(:,:,25), ncol, lchnk )
      call outfld( 'CPH26', cph_rate(:,:,26), ncol, lchnk )
      call outfld( 'CPH27', cph_rate(:,:,27), ncol, lchnk )
      call outfld( 'CPH28', cph_rate(:,:,28), ncol, lchnk )
      call outfld( 'CPH29', cph_rate(:,:,29), ncol, lchnk )
      call outfld( 'CPH30', cph_rate(:,:,30), ncol, lchnk )
      call outfld( 'CPH31', cph_rate(:,:,31), ncol, lchnk )
      call outfld( 'CPH32', cph_rate(:,:,32), ncol, lchnk )
      call outfld( 'CPH33', cph_rate(:,:,33), ncol, lchnk )
      call outfld( 'CPH34', cph_rate(:,:,34), ncol, lchnk )
      call outfld( 'CPH35', cph_rate(:,:,35), ncol, lchnk )
      call outfld( 'CPH36', cph_rate(:,:,36), ncol, lchnk )
      call outfld( 'CPH37', cph_rate(:,:,37), ncol, lchnk )
      call outfld( 'CPH38', cph_rate(:,:,38), ncol, lchnk )
      call outfld( 'CPH39', cph_rate(:,:,39), ncol, lchnk )
      call outfld( 'CPH40', cph_rate(:,:,40), ncol, lchnk )
      call outfld( 'CPH41', cph_rate(:,:,41), ncol, lchnk )
      call outfld( 'QCP', cph_tot(:,:), ncol, lchnk )

      end subroutine cph

      end module mo_cph
