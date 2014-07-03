
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
      use abortutils,  only : endrun
      use mo_tracname, only : solsym
      use chem_mods,   only : frc_from_dataset
      use shr_kind_mod,only : r8 => shr_kind_r8
      use cam_logfile, only : iulog

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios

      clscnt(:) = (/      1,     0,     0,    30,     0 /)

      cls_rxt_cnt(:,1) = (/      4,     1,     0,     1 /)
      cls_rxt_cnt(:,4) = (/      0,     9,    19,    30 /)

      solsym(: 31) = (/ 'O3              ','OH              ','HO2             ','H2O2            ','NO              ', &
                        'NO2             ','HNO3            ','CO              ','CH2O            ','CH3O2           ', &
                        'CH3OOH          ','DMS             ','SO2             ','ISOP            ','H2SO4           ', &
                        'SOAG            ','so4_a1          ','pom_a1          ','soa_a1          ','bc_a1           ', &
                        'dst_a1          ','ncl_a1          ','num_a1          ','so4_a2          ','soa_a2          ', &
                        'ncl_a2          ','num_a2          ','dst_a3          ','ncl_a3          ','so4_a3          ', &
                        'num_a3          ' /)

      adv_mass(: 31) = (/    47.998200_r8,    17.006800_r8,    33.006200_r8,    34.013600_r8,    30.006140_r8, &
                             46.005540_r8,    63.012340_r8,    28.010400_r8,    30.025200_r8,    47.032000_r8, &
                             48.039400_r8,    62.132400_r8,    64.064800_r8,    68.114200_r8,    98.078400_r8, &
                             12.011000_r8,   115.107340_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                            135.064039_r8,    58.442468_r8,     1.007400_r8,   115.107340_r8,    12.011000_r8, &
                             58.442468_r8,     1.007400_r8,   135.064039_r8,    58.442468_r8,   115.107340_r8, &
                              1.007400_r8 /)

      crb_mass(: 31) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                             12.011000_r8,    24.022000_r8,     0.000000_r8,    60.055000_r8,     0.000000_r8, &
                             12.011000_r8,     0.000000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8 /)

      fix_mass(:  5) = (/ 0.00000000_r8, 28.0134800_r8, 31.9988000_r8, 18.0142000_r8, 16.0406000_r8 /)

      clsmap(:  1,1) = (/    8 /)
      clsmap(: 30,4) = (/    1,   2,   3,   4,   5,   6,   7,   9,  10,  11, &
                            12,  13,  14,  15,  16,  17,  18,  19,  20,  21, &
                            22,  23,  24,  25,  26,  27,  28,  29,  30,  31 /)

      permute(: 30,4) = (/   30,  26,  27,  21,  29,  25,   1,  23,  28,  24, &
                             20,  19,  22,   2,   3,   4,   5,   6,   7,   8, &
                              9,  10,  11,  12,  13,  14,  15,  16,  17,  18 /)

      diag_map(: 30) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  18,  20,  23, &
                            26,  29,  35,  39,  44,  57,  66,  75,  83,  92 /)

      extfrc_lst(:  9) = (/ 'NO              ','NO2             ','SO2             ','so4_a1          ','so4_a2          ', &
                            'pom_a1          ','bc_a1           ','num_a1          ','num_a2          ' /)

      frc_from_dataset(:  9) = (/ .false., .true., .true., .true., .true., &
                                  .true., .true., .true., .true. /)

      inv_lst(:  5) = (/ 'M               ', 'N2              ', 'O2              ', 'H2O             ', 'CH4             ' /)

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
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'j2oh            ', 'jh2o2           ', 'jno2            ', 'jch2o_a         ', &
                                     'jch2o_b         ', 'jch3ooh         ', 'out6            ', 'usr_HO2_HO2     ', &
                                     'usr_oh_co       ', 'het_no2_h2o     ', 'usr_DMS_OH      ', 'usr_SO2_OH      ', &
                                     'isop_oh         ', 'isop_o3         ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   8,  10,  16,  23, &
                                       25,  26,  27,  28 /)
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
                              '                ', '                ' /)
      pht_alias_lst(:,2) = (/ 'jo3_a           ', '                ', '                ', '                ', &
                              '                ', '                ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
