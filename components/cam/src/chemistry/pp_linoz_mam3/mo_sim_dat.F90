
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

      clscnt(:) = (/      1,     0,     0,    20,     0 /)

      cls_rxt_cnt(:,1) = (/      0,     0,     0,     1 /)
      cls_rxt_cnt(:,4) = (/      1,     6,     0,    20 /)

      solsym(: 21) = (/ 'H2O2            ','H2SO4           ','SO2             ','DMS             ','SOAG            ', &
                        'so4_a1          ','pom_a1          ','soa_a1          ','bc_a1           ','dst_a1          ', &
                        'ncl_a1          ','num_a1          ','so4_a2          ','soa_a2          ','ncl_a2          ', &
                        'num_a2          ','dst_a3          ','ncl_a3          ','so4_a3          ','num_a3          ', &
                        'O3              ' /)

      adv_mass(: 21) = (/      34.013600_r8,      98.078400_r8,      64.064800_r8,      62.132400_r8,      12.011000_r8, &
                              115.107340_r8,      12.011000_r8,      12.011000_r8,      12.011000_r8,     135.064039_r8, &
                               58.442468_r8,       1.007400_r8,     115.107340_r8,      12.011000_r8,      58.442468_r8, &
                                1.007400_r8,     135.064039_r8,      58.442468_r8,     115.107340_r8,       1.007400_r8, &
                               47.998200_r8 /)

      crb_mass(: 21) = (/       0.000000_r8,       0.000000_r8,       0.000000_r8,      24.022000_r8,      12.011000_r8, &
                                0.000000_r8,      12.011000_r8,      12.011000_r8,      12.011000_r8,       0.000000_r8, &
                                0.000000_r8,       0.000000_r8,       0.000000_r8,      12.011000_r8,       0.000000_r8, &
                                0.000000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8,       0.000000_r8, &
                                0.000000_r8 /)

      fix_mass(:  8) = (/ 0.00000000_r8, 28.0134800_r8, 31.9988000_r8, 18.0142000_r8, 17.0068000_r8, &
                          62.0049400_r8, 33.0062000_r8, 47.9982000_r8 /)

      clsmap(:  1,1) = (/   21 /)
      clsmap(: 20,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  18,  19,  20 /)

      permute(: 20,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                             11,  12,  13,  14,  15,  16,  17,  18,  19,  20 /)

      diag_map(: 20) = (/    1,   2,   4,   6,   7,   8,   9,  10,  11,  12, &
                            13,  14,  15,  16,  17,  18,  19,  20,  21,  22 /)

      extfrc_lst(:  7) = (/ 'SO2             ','so4_a1          ','so4_a2          ','pom_a1          ','bc_a1           ', &
                            'num_a1          ','num_a2          ' /)

      frc_from_dataset(:  7) = (/ .true., .true., .true., .true., .true., &
                                  .true., .true. /)

      inv_lst(:  8) = (/ 'M               ', 'N2              ', 'O2              ', 'H2O             ', 'OH              ', &
                         'NO3             ', 'HO2             ', 'cnst_O3         ' /)

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
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jh2o2           ', 'usr_HO2_HO2     ', 'usr_SO2_OH      ', 'usr_DMS_OH      ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   4,   6 /)
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
      pht_alias_lst(:,1) = (/ '                ' /)
      pht_alias_lst(:,2) = (/ '                ' /)
      pht_alias_mult(:,1) = (/ 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
