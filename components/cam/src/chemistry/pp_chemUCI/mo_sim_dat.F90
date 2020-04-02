
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

      clscnt(:) = (/      7,     0,     0,    24,     0 /)

      cls_rxt_cnt(:,1) = (/      9,    10,     0,     7 /)
      cls_rxt_cnt(:,4) = (/      2,    30,    49,    24 /)

      solsym(: 31) = (/ 'O3              ','OH              ','HO2             ','H2O2            ','CH2O            ', &
                        'CH3O2           ','CH3OOH          ','NO              ','NO2             ','NO3             ', &
                        'N2O5            ','HNO3            ','HNO4            ','PAN             ','CO              ', &
                        'CH4             ','C2H6            ','C3H8            ','C2H4            ','ROHO2           ', &
                        'CH3COCH3        ','C2H5O2          ','C2H5OOH         ','CH3CHO          ','CH3CO3          ', &
                        'ISOP            ','ISOPO2          ','MVKMACR         ','MVKO2           ','e90             ', &
                        'O3S             ' /)

      adv_mass(: 31) = (/      47.998200_r8,      17.006800_r8,      33.006200_r8,      34.013600_r8,      30.025200_r8, &
                               47.032000_r8,      48.039400_r8,      30.006140_r8,      46.005540_r8,      62.004940_r8, &
                              108.010480_r8,      63.012340_r8,      79.011740_r8,     121.047940_r8,      28.010400_r8, &
                               16.040600_r8,      30.066400_r8,      44.092200_r8,      28.051600_r8,      77.057200_r8, &
                               58.076800_r8,      61.057800_r8,      62.065200_r8,      44.051000_r8,      75.042400_r8, &
                               68.114200_r8,     117.119800_r8,      70.087800_r8,     119.093400_r8,      47.998200_r8, &
                               47.998200_r8 /)

      crb_mass(: 31) = (/       0.000000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8,      12.011000_r8, &
                               12.011000_r8,      12.011000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8, &
                                0.000000_r8,       0.000000_r8,       0.000000_r8,      24.022000_r8,      12.011000_r8, &
                               12.011000_r8,      24.022000_r8,      36.033000_r8,      24.022000_r8,      24.022000_r8, &
                               36.033000_r8,      24.022000_r8,      24.022000_r8,      24.022000_r8,      24.022000_r8, &
                               60.055000_r8,      60.055000_r8,      48.044000_r8,      48.044000_r8,       0.000000_r8, &
                                0.000000_r8 /)

      fix_mass(:  5) = (/ 0.00000000_r8, 28.0134800_r8, 31.9988000_r8, 18.0142000_r8, 2.01480000_r8 /)

      clsmap(:  7,1) = (/   15,  16,  17,  18,  21,  30,  31 /)
      clsmap(: 24,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  22,  23,  24,  25,  19,  20, &
                            26,  27,  28,  29 /)

      permute(: 24,4) = (/   19,  22,  24,   1,  10,  20,   6,  23,  18,  21, &
                              3,   7,   8,   2,  13,   9,  16,  17,   4,  12, &
                              5,  14,  15,  11 /)

      diag_map(: 24) = (/    1,   4,   7,  11,  17,  23,  28,  32,  37,  42, &
                            47,  56,  65,  75,  86,  96, 105, 116, 132, 148, &
                           160, 183, 199, 219 /)

      extfrc_lst(:  2) = (/ 'NO              ','NO2             ' /)

      frc_from_dataset(:  2) = (/ .false., .true. /)

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
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo1d            ', 'jo2             ', 'jh2o2           ', 'jch2o_a         ', &
                                     'jch2o_b         ', 'jch3ooh         ', 'jc2h5ooh        ', 'jno2            ', &
                                     'jno3            ', 'jn2o5           ', 'jhno3           ', 'jho2no2         ', &
                                     'jch3cho         ', 'jpan            ', 'jacet           ', 'jmvk            ', &
                                     'uci1            ', 'uci2            ', 'uci3            ', 'uci4            ', &
                                     'uci5            ', 'uci6            ', 'HNO4            ', 'N2O5            ', &
                                     'PAN             ', 'uci7            ', 'uci8            ', 'uci9            ', &
                                     'ucih1           ', 'ucih2           ', 'ucih3           ', 'usr_e90         ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                                       11,  12,  13,  14,  15,  16,  17,  18,  19,  35, &
                                       36,  45,  47,  48,  49,  50,  51,  52,  80,  81, &
                                       82,  83 /)
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
      pht_alias_lst(:,2) = (/ 'jo3_a           ', '                ', '                ', '                ', &
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
