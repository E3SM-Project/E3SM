
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

      clscnt(:) = (/     17,     0,     0,    49,     0 /)

      cls_rxt_cnt(:,1) = (/     12,    41,     0,    17 /)
      cls_rxt_cnt(:,4) = (/     23,   102,   129,    49 /)

      solsym(: 66) = (/ 'O3              ','O               ','O1D             ','O2              ','O2_1S           ', &
                        'O2_1D           ','N2O             ','N               ','NO              ','NO2             ', &
                        'NO3             ','HNO3            ','HO2NO2          ','N2O5            ','CH4             ', &
                        'CH3O2           ','CH3OOH          ','CH2O            ','CO              ','H2              ', &
                        'H               ','OH              ','HO2             ','H2O2            ','CLY             ', &
                        'BRY             ','CL              ','CL2             ','CLO             ','OCLO            ', &
                        'CL2O2           ','HCL             ','HOCL            ','CLONO2          ','BRCL            ', &
                        'BR              ','BRO             ','HBR             ','HOBR            ','BRONO2          ', &
                        'CH3CL           ','CH3BR           ','CFC11           ','CFC12           ','CFC113          ', &
                        'HCFC22          ','CCL4            ','CH3CCL3         ','CF3BR           ','CF2CLBR         ', &
                        'CO2             ','N2p             ','O2p             ','Np              ','Op              ', &
                        'NOp             ','e               ','N2D             ','OCS             ','S               ', &
                        'SO              ','SO2             ','SO3             ','HSO3            ','H2SO4           ', &
                        'H2O             ' /)

      adv_mass(: 66) = (/    47.998200_r8,    15.999400_r8,    15.999400_r8,    31.998800_r8,    31.998800_r8, &
                             31.998800_r8,    44.012880_r8,    14.006740_r8,    30.006140_r8,    46.005540_r8, &
                             62.004940_r8,    63.012340_r8,    79.011740_r8,   108.010480_r8,    16.040600_r8, &
                             47.032000_r8,    48.039400_r8,    30.025200_r8,    28.010400_r8,     2.014800_r8, &
                              1.007400_r8,    17.006800_r8,    33.006200_r8,    34.013600_r8,   100.916850_r8, &
                             99.716850_r8,    35.452700_r8,    70.905400_r8,    51.452100_r8,    67.451500_r8, &
                            102.904200_r8,    36.460100_r8,    52.459500_r8,    97.457640_r8,   115.356700_r8, &
                             79.904000_r8,    95.903400_r8,    80.911400_r8,    96.910800_r8,   141.908940_r8, &
                             50.485900_r8,    94.937200_r8,   137.367503_r8,   120.913206_r8,   187.375310_r8, &
                             86.467906_r8,   153.821800_r8,   133.402300_r8,   148.910210_r8,   165.364506_r8, &
                             44.009800_r8,    28.013480_r8,    31.998800_r8,    14.006740_r8,    15.999400_r8, &
                             30.006140_r8, 0.548567E-03_r8,    14.006740_r8,    60.076400_r8,    32.066000_r8, &
                             48.065400_r8,    64.064800_r8,    80.064200_r8,    81.071600_r8,    98.078400_r8, &
                             18.014200_r8 /)

      crb_mass(: 66) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    24.022000_r8, &
                             12.011000_r8,    12.011000_r8,    24.022000_r8,    12.011000_r8,    12.011000_r8, &
                             12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8 /)

      fix_mass(:  2) = (/ 0.00000000_r8, 28.0134800_r8 /)

      clsmap(: 17,1) = (/   15,   7,  19,  20,  41,  42,  43,  44,  45,  46, &
                            47,  48,  49,  50,  51,  25,  26 /)
      clsmap(: 49,4) = (/    1,   2,   3,   4,   5,   6,   8,   9,  10,  22, &
                            11,  12,  13,  14,  16,  17,  18,  21,  23,  24, &
                            66,  27,  28,  29,  30,  31,  32,  33,  34,  35, &
                            36,  37,  38,  39,  40,  52,  53,  54,  55,  56, &
                            58,  57,  59,  60,  61,  62,  63,  64,  65 /)

      permute(: 49,4) = (/   42,  39,  45,  49,   3,   2,  27,  38,  40,  43, &
                             33,  22,  13,   8,  28,  11,  46,  32,  47,  18, &
                             44,  48,   5,  37,   9,   1,  36,  29,  30,   6, &
                             35,  41,  20,  23,  17,  19,  24,  14,  15,  25, &
                             16,  26,  10,  21,  34,  31,   7,  12,   4 /)

      diag_map(: 49) = (/    1,   4,   7,   9,  12,  14,  18,  22,  28,  33, &
                            40,  47,  53,  60,  67,  73,  78,  86,  95, 104, &
                           110, 117, 123, 131, 140, 151, 162, 169, 179, 189, &
                           199, 206, 215, 228, 242, 261, 286, 308, 347, 373, &
                           397, 419, 450, 470, 490, 510, 535, 557, 588 /)

      extfrc_lst(: 12) = (/ 'NO              ','NO2             ','CO              ','SO2             ','Op              ', &
                            'O2p             ','Np              ','N2p             ','N2D             ','N               ', &
                            'e               ','OH              ' /)

      frc_from_dataset(: 12) = (/ .true., .true., .true., .true., .false., &
                                  .false., .false., .false., .false., .false., &
                                  .false., .false. /)

      inv_lst(:  2) = (/ 'M               ', 'N2              ' /)

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
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo2_a           ', 'jo2_b           ', 'jo3_a           ', 'jo3_b           ', &
                                     'jn2o            ', 'jno             ', 'jno_i           ', 'jno2            ', &
                                     'jn2o5_a         ', 'jn2o5_b         ', 'jhno3           ', 'jno3_a          ', &
                                     'jno3_b          ', 'jho2no2_a       ', 'jho2no2_b       ', 'jch3ooh         ', &
                                     'jch2o_a         ', 'jch2o_b         ', 'jh2o_a          ', 'jh2o_b          ', &
                                     'jh2o_c          ', 'jh2o2           ', 'jcl2            ', 'jclo            ', &
                                     'joclo           ', 'jcl2o2          ', 'jhocl           ', 'jhcl            ', &
                                     'jclono2_a       ', 'jclono2_b       ', 'jbrcl           ', 'jbro            ', &
                                     'jhobr           ', 'jbrono2_a       ', 'jbrono2_b       ', 'jch3cl          ', &
                                     'jccl4           ', 'jch3ccl3        ', 'jcfcl3          ', 'jcf2cl2         ', &
                                     'jcfc113         ', 'jhcfc22         ', 'jch3br          ', 'jcf3br          ', &
                                     'jcf2clbr        ', 'jco2            ', 'jch4_a          ', 'jch4_b          ', &
                                     'jh2so4          ', 'jso2            ', 'jso3            ', 'jocs            ', &
                                     'jso             ', 'jeuv_1          ', 'jeuv_2          ', 'jeuv_3          ', &
                                     'jeuv_4          ', 'jeuv_5          ', 'jeuv_6          ', 'jeuv_7          ', &
                                     'jeuv_8          ', 'jeuv_9          ', 'jeuv_10         ', 'jeuv_11         ', &
                                     'jeuv_12         ', 'jeuv_13         ', 'jeuv_14         ', 'jeuv_15         ', &
                                     'jeuv_16         ', 'jeuv_17         ', 'jeuv_18         ', 'jeuv_19         ', &
                                     'jeuv_20         ', 'jeuv_21         ', 'jeuv_22         ', 'jeuv_23         ', &
                                     'jeuv_24         ', 'jeuv_25         ', 'jeuv_26         ', 'usr_O_O2        ', &
                                     'cph_O_O3        ', 'usr_O_O         ', 'cph_O2_1S_O     ', 'cph_O2_1S_O2    ', &
                                     'cph_O2_1S_N2    ', 'cph_O2_1S_O3    ', 'ag2             ', 'cph_O2_1D_O     ', &
                                     'cph_O2_1D_O2    ', 'cph_O2_1D_N2    ', 'ag1             ', 'cph_O1D_N2      ', &
                                     'cph_O1D_O2      ', 'cph_O1D_O2b     ', 'cph_N2D_O2      ', 'cph_N2D_O       ', &
                                     'cph_N_O2        ', 'cph_N_NO        ', 'cph_NO_HO2      ', 'cph_NO_O3       ', &
                                     'cph_NO2_O       ', 'tag_NO2_NO3     ', 'usr_N2O5_M      ', 'tag_NO2_OH      ', &
                                     'usr_HNO3_OH     ', 'tag_NO2_HO2     ', 'usr_HO2NO2_M    ', 'usr_CO_OH_b     ', &
                                     'cph_H_O2        ', 'cph_H_O3        ', 'cph_H_HO2       ', 'cph_OH_O        ', &
                                     'cph_OH_O3       ', 'cph_OH_HO2      ', 'cph_HO2_O       ', 'cph_HO2_O3      ', &
                                     'usr_HO2_HO2     ', 'tag_CLO_CLO     ', 'usr_CL2O2_M     ', 'usr_SO3_H2O     ', &
                                     'het1            ', 'het2            ', 'het3            ', 'het4            ', &
                                     'het5            ', 'het6            ', 'het7            ', 'het8            ', &
                                     'het9            ', 'het10           ', 'het11           ', 'het12           ', &
                                     'het13           ', 'het14           ', 'het15           ', 'het16           ', &
                                     'het17           ', 'ion_Op_O2       ', 'ion_Op_N2       ', 'ion_N2p_Oa      ', &
                                     'ion_O2p_N       ', 'ion_O2p_NO      ', 'ion_Np_O2a      ', 'ion_Np_O2b      ', &
                                     'ion_Np_O        ', 'ion_N2p_O2      ', 'ion_N2p_Ob      ', 'elec1           ', &
                                     'elec2           ', 'elec3           ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                                       11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                                       21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                                       31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                                       41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                                       51,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                                       61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                                       71,  72,  73,  74,  75,  76,  77,  78,  79,  80, &
                                       81,  82,  83,  84,  85,  86,  88,  89,  90,  91, &
                                       92,  93,  94,  95, 114, 115, 116, 117, 120, 121, &
                                      122, 125, 126, 127, 128, 133, 135, 144, 145, 146, &
                                      148, 150, 151, 152, 157, 158, 159, 177, 178, 222, &
                                      223, 224, 225, 226, 227, 228, 229, 230, 231, 232, &
                                      233, 234, 235, 236, 237, 238, 239, 240, 241, 242, &
                                      244, 245, 246, 247, 248, 249, 251, 252, 253, 254 /)
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
      pht_alias_lst(:,1) = (/ 'userdefined     ', 'userdefined     ', '                ', '                ', &
                              '                ', 'userdefined     ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ' /)
      pht_alias_lst(:,2) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
