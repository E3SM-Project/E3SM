
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

      clscnt(:) = (/     10,     0,     0,    43,     0 /)

      cls_rxt_cnt(:,1) = (/      9,    12,     0,    10 /)
      cls_rxt_cnt(:,4) = (/      3,    30,    53,    43 /)

      solsym(: 53) = (/ 'O3              ','OH              ','HO2             ','H2O2            ','CH2O            ', &
                        'CH3O2           ','CH3OOH          ','NO              ','NO2             ','NO3             ', &
                        'N2O5            ','HNO3            ','HNO4            ','PAN             ','CO              ', &
                        'CH4             ','C2H6            ','C3H8            ','C2H4            ','CH3COCH3        ', &
                        'C2H5O2          ','C2H5OOH         ','CH3CHO          ','CH3CO3          ','ISOP            ', &
                        'ISOPO2          ','MVKMACR         ','MVKO2           ','ROHO2           ','e90             ', &
                        'O3S             ','N2O             ','NOYS            ','ASAD            ','SO2             ', &
                        'DMS             ','H2SO4           ','SOAG            ','so4_a1          ','pom_a1          ', &
                        'soa_a1          ','bc_a1           ','dst_a1          ','ncl_a1          ','num_a1          ', &
                        'so4_a2          ','soa_a2          ','ncl_a2          ','num_a2          ','dst_a3          ', &
                        'ncl_a3          ','so4_a3          ','num_a3          ' /)

      adv_mass(: 53) = (/      47.998200_r8,      17.006800_r8,      33.006200_r8,      34.013600_r8,      30.025200_r8, &
                               47.032000_r8,      48.039400_r8,      30.006140_r8,      46.005540_r8,      62.004940_r8, &
                              108.010480_r8,      63.012340_r8,      79.011740_r8,     121.047940_r8,      28.010400_r8, &
                               16.040600_r8,      30.066400_r8,      44.092200_r8,      28.051600_r8,      58.076800_r8, &
                               61.057800_r8,      62.065200_r8,      44.051000_r8,      75.042400_r8,      68.114200_r8, &
                              117.119800_r8,     102.050800_r8,     122.038600_r8,      49.005600_r8,   0.493710E-01_r8, &
                               80.064200_r8,      44.012880_r8,     150.977990_r8,      32.066000_r8,      64.064800_r8, &
                               62.132400_r8,      98.078400_r8,      12.011000_r8,     115.107340_r8,      12.011000_r8, &
                               12.011000_r8,      12.011000_r8,     135.064039_r8,      58.442468_r8,       1.007400_r8, &
                              115.107340_r8,      12.011000_r8,      58.442468_r8,       1.007400_r8,     135.064039_r8, &
                               58.442468_r8,     115.107340_r8,       1.007400_r8 /)

      crb_mass(: 53) = (/       0.000000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8,      12.011000_r8, &
                               12.011000_r8,      12.011000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8, &
                                0.000000_r8,       0.000000_r8,       0.000000_r8,      24.022000_r8,      12.011000_r8, &
                               12.011000_r8,      24.022000_r8,      36.033000_r8,      24.022000_r8,      36.033000_r8, &
                               24.022000_r8,      24.022000_r8,      24.022000_r8,      24.022000_r8,      60.055000_r8, &
                               60.055000_r8,      12.011000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8, &
                                0.000000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8, &
                               24.022000_r8,       0.000000_r8,      12.011000_r8,       0.000000_r8,      12.011000_r8, &
                               12.011000_r8,      12.011000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8, &
                                0.000000_r8,      12.011000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8, &
                                0.000000_r8,       0.000000_r8,       0.000000_r8 /)

      fix_mass(:  5) = (/ 0.00000000_r8, 28.0134800_r8, 31.9988000_r8, 18.0142000_r8, 2.01480000_r8 /)

      clsmap(: 10,1) = (/   15,  16,  17,  18,  20,  30,  31,  32,  33,  34 /)
      clsmap(: 43,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  21,  22,  23,  24,  19,  29, &
                            25,  26,  27,  28,  36,  35,  37,  38,  39,  40, &
                            41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                            51,  52,  53 /)

      permute(: 43,4) = (/   39,  42,  40,  20,  29,  41,  25,  37,  43,  38, &
                             22,  26,  27,  21,  32,  28,  35,  36,  23,  31, &
                             24,  33,  34,  30,  19,  18,   1,   2,   3,   4, &
                              5,   6,   7,   8,   9,  10,  11,  12,  13,  14, &
                             15,  16,  17 /)

      diag_map(: 43) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  19,  22,  25, &
                            28,  31,  35,  41,  47,  52,  56,  61,  66,  71, &
                            80,  89,  99, 110, 120, 129, 143, 155, 171, 191, &
                           207, 233, 246 /)

      extfrc_lst(: 10) = (/ 'NO              ','NO2             ','SO2             ','so4_a1          ','so4_a2          ', &
                            'pom_a1          ','bc_a1           ','num_a1          ','num_a2          ','ASAD            ' /)

      frc_from_dataset(: 10) = (/ .false., .true., .true., .true., .true., &
                                  .true., .true., .true., .true., .true. /)

      inv_lst(:  5) = (/ 'M               ', 'N2              ', 'O2              ', 'H2O             ', 'H2              ' /)

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
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo2             ', 'jh2o2           ', 'jch2o_a         ', 'jch2o_b         ', &
                                     'jch3ooh         ', 'jc2h5ooh        ', 'jno2            ', 'jno3            ', &
                                     'jn2o5           ', 'jhno3           ', 'jhno4           ', 'jch3cho         ', &
                                     'jpan            ', 'jacet_a         ', 'jacet_b         ', 'jmvkmacr        ', &
                                     'usr_jo3_a       ', 'usr_jo3_b       ', 'usr_jo3_c       ', '3b_H2O2         ', &
                                     '2b_H2O2         ', '3b_HNO3         ', '3b_HNO4         ', '3b_N2O5         ', &
                                     '3b_PAN          ', 'Keq_HNO4        ', 'Keq_N2O5        ', 'Keq_PAN         ', &
                                     'het_N2O5        ', 'het_NO3         ', 'het_HO2         ', '3b_OHHNO3       ', &
                                     '3b_OH_CO        ', 'usr_DMS_OH      ', 'usr_SO2_OH      ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                                       11,  12,  13,  14,  15,  16,  17,  18,  19,  24, &
                                       25,  33,  34,  35,  36,  37,  38,  39,  40,  41, &
                                       42,  44,  48,  85,  86 /)
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
      pht_alias_lst(:,1) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ' /)
      pht_alias_lst(:,2) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
