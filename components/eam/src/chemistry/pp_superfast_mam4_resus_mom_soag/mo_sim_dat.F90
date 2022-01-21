
      module mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat

      use chem_mods,   only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass, crb_mass
      use chem_mods,   only : diag_map
      use chem_mods,   only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods,   only : pht_alias_lst, pht_alias_mult
      use chem_mods,   only : extfrc_lst, inv_lst, slvd_lst
      use cam_abortutils,  only : endrun
      use mo_tracname, only : solsym
      use chem_mods,   only : frc_from_dataset
      use shr_kind_mod,only : r8 => shr_kind_r8
      use cam_logfile, only : iulog

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios

      clscnt(:) = (/      6,     0,     0,    53,     0 /)

      cls_rxt_cnt(:,1) = (/      9,     8,     0,     6 /)
      cls_rxt_cnt(:,4) = (/      2,    30,    54,    53 /)

      solsym(: 59) = (/ 'O3              ','OH              ','HO2             ','H2O2            ','CH2O            ', &
                        'CH3O2           ','CH3OOH          ','NO              ','NO2             ','NO3             ', &
                        'N2O5            ','HNO3            ','HO2NO2          ','PAN             ','CO              ', &
                        'C2H6            ','C3H8            ','C2H4            ','ROHO2           ','CH3COCH3        ', &
                        'C2H5O2          ','C2H5OOH         ','CH3CHO          ','CH3CO3          ','ISOP            ', &
                        'ISOPO2          ','MVKMACR         ','MVKO2           ','E90             ','O3LNZ           ', &
                        'DMS             ','SO2             ','H2SO4           ','SOAG            ','so4_a1          ', &
                        'so4_a2          ','so4_a3          ','pom_a1          ','pom_a3          ','pom_a4          ', &
                        'soa_a1          ','soa_a2          ','soa_a3          ','bc_a1           ','bc_a3           ', &
                        'bc_a4           ','dst_a1          ','dst_a3          ','ncl_a1          ','ncl_a2          ', &
                        'ncl_a3          ','mom_a1          ','mom_a2          ','mom_a3          ','mom_a4          ', &
                        'num_a1          ','num_a2          ','num_a3          ','num_a4          ' /)

      adv_mass(: 59) = (/      47.998200_r8,      17.006800_r8,      33.006200_r8,      34.013600_r8,      30.025200_r8, &
                               47.032000_r8,      48.039400_r8,      30.006140_r8,      46.005540_r8,      62.004940_r8, &
                              108.010480_r8,      63.012340_r8,      79.011740_r8,     121.047940_r8,      28.010400_r8, &
                               30.066400_r8,      44.092200_r8,      28.051600_r8,      77.057200_r8,      58.076800_r8, &
                               61.057800_r8,      62.065200_r8,      44.051000_r8,      75.042400_r8,      68.114200_r8, &
                              117.119800_r8,      70.087800_r8,     119.093400_r8,      47.998200_r8,      47.998200_r8, &
                               62.132400_r8,      64.064800_r8,      98.078400_r8,      12.011000_r8,     115.107340_r8, &
                              115.107340_r8,     115.107340_r8,      12.011000_r8,      12.011000_r8,      12.011000_r8, &
                               12.011000_r8,      12.011000_r8,      12.011000_r8,      12.011000_r8,      12.011000_r8, &
                               12.011000_r8,     135.064039_r8,     135.064039_r8,      58.442468_r8,      58.442468_r8, &
                               58.442468_r8,  250092.672000_r8,  250092.672000_r8,  250092.672000_r8,  250092.672000_r8, &
                                1.007400_r8,       1.007400_r8,       1.007400_r8,       1.007400_r8 /)

      crb_mass(: 59) = (/       0.000000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8,      12.011000_r8, &
                               12.011000_r8,      12.011000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8, &
                                0.000000_r8,       0.000000_r8,       0.000000_r8,      24.022000_r8,      12.011000_r8, &
                               24.022000_r8,      36.033000_r8,      24.022000_r8,      24.022000_r8,      36.033000_r8, &
                               24.022000_r8,      24.022000_r8,      24.022000_r8,      24.022000_r8,      60.055000_r8, &
                               60.055000_r8,      48.044000_r8,      48.044000_r8,       0.000000_r8,       0.000000_r8, &
                               24.022000_r8,       0.000000_r8,       0.000000_r8,      12.011000_r8,       0.000000_r8, &
                                0.000000_r8,       0.000000_r8,      12.011000_r8,      12.011000_r8,      12.011000_r8, &
                               12.011000_r8,      12.011000_r8,      12.011000_r8,      12.011000_r8,      12.011000_r8, &
                               12.011000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8, &
                                0.000000_r8,  102333.720000_r8,  102333.720000_r8,  102333.720000_r8,  102333.720000_r8, &
                                0.000000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8 /)

      fix_mass(:  6) = (/ 0.00000000_r8, 28.0134800_r8, 31.9988000_r8, 18.0142000_r8, 2.01480000_r8, &
                          16.0406000_r8 /)

      clsmap(:  6,1) = (/   15,  16,  17,  20,  29,  30 /)
      clsmap(: 53,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  21,  22,  23,  24,  18,  19, &
                            25,  26,  27,  28,  31,  32,  33,  34,  35,  36, &
                            37,  38,  39,  40,  41,  42,  43,  44,  45,  46, &
                            47,  48,  49,  50,  51,  52,  53,  54,  55,  56, &
                            57,  58,  59 /)

      permute(: 53,4) = (/   53,  52,  50,  29,  39,  51,  37,  47,  49,  48, &
                             31,  35,  36,  30,  42,  38,  45,  46,  33,  41, &
                             34,  43,  44,  40,  32,  28,   1,   2,   3,   4, &
                              5,   6,   7,   8,   9,  10,  11,  12,  13,  14, &
                             15,  16,  17,  18,  19,  20,  21,  22,  23,  24, &
                             25,  26,  27 /)

      diag_map(: 53) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                            21,  22,  23,  24,  25,  26,  27,  29,  31,  34, &
                            37,  43,  48,  54,  60,  64,  69,  74,  79,  84, &
                            93, 102, 112, 123, 133, 142, 156, 170, 183, 203, &
                           219, 245, 261 /)

      extfrc_lst(: 11) = (/ 'NO              ','NO2             ','SO2             ','so4_a1          ','so4_a2          ', &
                            'pom_a4          ','bc_a4           ','num_a1          ','num_a2          ','num_a4          ', &
                            'SOAG            ' /)

      frc_from_dataset(: 11) = (/ .false., .true., .true., .true., .true., &
                                  .true., .true., .true., .true., .true., &
                                  .true. /)

      inv_lst(:  6) = (/ 'M               ', 'N2              ', 'O2              ', 'H2O             ', 'H2              ', &
                         'CH4             ' /)

      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo1dU           ', 'jo2_b           ', 'jh2o2           ', 'jch2o_a         ', &
                                     'jch2o_b         ', 'jch3ooh         ', 'jc2h5ooh        ', 'jno2            ', &
                                     'jno3_a          ', 'jno3_b          ', 'jn2o5_a         ', 'jn2o5_b         ', &
                                     'jhno3           ', 'jho2no2_a       ', 'jho2no2_b       ', 'jch3cho         ', &
                                     'jpan            ', 'jacet           ', 'jmvk            ', 'uci1            ', &
                                     'uci2            ', 'uci3            ', 'lco_h           ', 'lco_ho2         ', &
                                     'lh2_ho2         ', 'lch4            ', 'lc2h6           ', 'lc3h8           ', &
                                     'lc2h4_oh        ', 'lc2h4_o3        ', 'lisop_o3        ', 'lisop_oh        ', &
                                     'lch2o           ', 'lo3_oh          ', 'po3_oh          ', 'lo3_ho2         ', &
                                     'lho2_oh         ', 'uci4            ', 'uci5            ', 'ph2o2           ', &
                                     'lh2o2           ', 'lo3_no          ', 'lno_ho2         ', 'lo3_no2         ', &
                                     'lno3_oh         ', 'lno3_no         ', 'lhno4           ', 'lhno3           ', &
                                     'uci6            ', 'lno2_oh         ', 'HO2NO2f         ', 'N2O5f           ', &
                                     'PANf            ', 'uci7            ', 'uci8            ', 'uci9            ', &
                                     'lch3o2_ho2      ', 'lch3o2_no       ', 'lch3o2          ', 'lch3ooh         ', &
                                     'lc2h5o2_no      ', 'lc2h5o2         ', 'lc2h5o2_ch3     ', 'lc2h5o2_ho2     ', &
                                     'lc2h5ooh_a      ', 'lc2h5ooh_b      ', 'lch3cho_oh      ', 'lch3cho_no3     ', &
                                     'lch3co3_no      ', 'lch3co3_ch3     ', 'lch3co3         ', 'lch3coch3_a     ', &
                                     'lch3coch3_b     ', 'lroho2_no       ', 'lroho2_ho2      ', 'lroho2_ch3o2    ', &
                                     'lisopo2_no      ', 'lisopo2_ho2     ', 'lisopo2_ch3     ', 'lmvkmacr_o3     ', &
                                     'lmvkmacr_oh     ', 'lmvko2_no       ', 'lmvko2_ho2      ', 'usr_e90         ', &
                                     'ldms_oh         ', 'usr_DMS_OH      ', 'usr_SO2_OH      ', 'ldms_no3        ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                                       11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                                       21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                                       31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                                       41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                                       51,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                                       61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                                       71,  72,  73,  74,  75,  76,  77,  78,  79,  80, &
                                       81,  82,  83,  84,  85,  86,  87,  88 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_lst; error = ',ios
         call endrun
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_mult; error = ',ios
         call endrun
      end if
      pht_alias_lst(:,1) = (/ '                ', 'userdefined     ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ' /)
      pht_alias_lst(:,2) = (/ 'jo3_a           ', '                ', '                ', '                ', &
                              '                ', '                ', 'jch3ooh         ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
