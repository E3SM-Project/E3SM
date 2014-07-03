
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

      clscnt(:) = (/      0,     0,     0,    16,     0 /)

      cls_rxt_cnt(:,4) = (/      0,     2,     0,    16 /)

      solsym(: 16) = (/ 'H2O2            ','SO2             ','SO4             ','DMS             ','OC1             ', &
                        'OC2             ','CB1             ','CB2             ','SSLT01          ','SSLT02          ', &
                        'SSLT03          ','SSLT04          ','DST01           ','DST02           ','DST03           ', &
                        'DST04           ' /)

      adv_mass(: 16) = (/    34.013600_r8,    64.064800_r8,    96.063600_r8,    62.132400_r8,    12.011000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    58.442468_r8,    58.442468_r8, &
                             58.442468_r8,    58.442468_r8,   135.064039_r8,   135.064039_r8,   135.064039_r8, &
                            135.064039_r8 /)

      crb_mass(: 16) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8,    24.022000_r8,    12.011000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8 /)

      fix_mass(:  7) = (/ 0.00000000_r8, 28.0134800_r8, 31.9988000_r8, 47.9982000_r8, 17.0068000_r8, &
                          62.0049400_r8, 33.0062000_r8 /)

      clsmap(: 16,4) = (/    1,   2,   3,   4,   7,   8,   5,   6,   9,  10, &
                            11,  12,  13,  14,  15,  16 /)

      permute(: 16,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                             11,  12,  13,  14,  15,  16 /)

      diag_map(: 16) = (/    1,   2,   3,   4,   5,   7,   8,  10,  11,  12, &
                            13,  14,  15,  16,  17,  18 /)

      extfrc_lst(:  2) = (/ 'SO2             ','SO4             ' /)

      frc_from_dataset(:  2) = (/ .true., .true. /)

      inv_lst(:  7) = (/ 'M               ', 'N2              ', 'O2              ', 'O3              ', 'OH              ', &
                         'NO3             ', 'HO2             ' /)

      end subroutine set_sim_dat

      end module mo_sim_dat
