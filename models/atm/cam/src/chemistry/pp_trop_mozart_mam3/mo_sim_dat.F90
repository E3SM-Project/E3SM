
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

      clscnt(:) = (/      4,     0,     0,   104,     0 /)

      cls_rxt_cnt(:,1) = (/      1,     8,     0,     4 /)
      cls_rxt_cnt(:,4) = (/      2,    62,   164,   104 /)

      solsym(:108) = (/ 'O3              ','O               ','O1D             ','N2O             ','NO              ', &
                        'NO2             ','NO3             ','HNO3            ','HO2NO2          ','N2O5            ', &
                        'H2              ','OH              ','HO2             ','H2O2            ','H               ', &
                        'CH4             ','CO              ','CH3O2           ','CH3OOH          ','CH2O            ', &
                        'CH3OH           ','C2H5OH          ','C2H4            ','EO              ','EO2             ', &
                        'CH3COOH         ','GLYALD          ','EOOH            ','C2H6            ','C2H5O2          ', &
                        'C2H5OOH         ','CH3CHO          ','CH3CO3          ','CH3COOOH        ','C3H6            ', &
                        'C3H8            ','C3H7O2          ','C3H7OOH         ','PO2             ','POOH            ', &
                        'CH3COCH3        ','RO2             ','ROOH            ','BIGENE          ','ENEO2           ', &
                        'MEK             ','MEKO2           ','MEKOOH          ','BIGALK          ','ALKO2           ', &
                        'ALKOOH          ','ISOP            ','ISOPO2          ','ISOPOOH         ','MVK             ', &
                        'MACR            ','MACRO2          ','MACROOH         ','MCO3            ','HYDRALD         ', &
                        'HYAC            ','CH3COCHO        ','XO2             ','XOOH            ','C10H16          ', &
                        'TERPO2          ','TERPOOH         ','TOLUENE         ','CRESOL          ','TOLO2           ', &
                        'TOLOOH          ','XOH             ','BIGALD          ','GLYOXAL         ','PAN             ', &
                        'ONIT            ','MPAN            ','ISOPNO3         ','ONITR           ','SO2             ', &
                        'DMS             ','NH3             ','NH4             ','Rn              ','Pb              ', &
                        'HCN             ','CH3CN           ','C2H2            ','HCOOH           ','HOCH2OO         ', &
                        'H2SO4           ','SOAG            ','SOA             ','so4_a1          ','pom_a1          ', &
                        'soa_a1          ','bc_a1           ','dst_a1          ','ncl_a1          ','num_a1          ', &
                        'so4_a2          ','soa_a2          ','ncl_a2          ','num_a2          ','dst_a3          ', &
                        'ncl_a3          ','so4_a3          ','num_a3          ' /)

      adv_mass(:108) = (/    47.998200_r8,    15.999400_r8,    15.999400_r8,    44.012880_r8,    30.006140_r8, &
                             46.005540_r8,    62.004940_r8,    63.012340_r8,    79.011740_r8,   108.010480_r8, &
                              2.014800_r8,    17.006800_r8,    33.006200_r8,    34.013600_r8,     1.007400_r8, &
                             16.040600_r8,    28.010400_r8,    47.032000_r8,    48.039400_r8,    30.025200_r8, &
                             32.040000_r8,    46.065800_r8,    28.051600_r8,    61.057800_r8,    77.057200_r8, &
                             60.050400_r8,    60.050400_r8,    78.064600_r8,    30.066400_r8,    61.057800_r8, &
                             62.065200_r8,    44.051000_r8,    75.042400_r8,    76.049800_r8,    42.077400_r8, &
                             44.092200_r8,    75.083600_r8,    76.091000_r8,    91.083000_r8,    92.090400_r8, &
                             58.076800_r8,    89.068200_r8,    90.075600_r8,    56.103200_r8,   105.108800_r8, &
                             72.102600_r8,   103.094000_r8,   104.101400_r8,    72.143800_r8,   103.135200_r8, &
                            104.142600_r8,    68.114200_r8,   117.119800_r8,   118.127200_r8,    70.087800_r8, &
                             70.087800_r8,   119.093400_r8,   120.100800_r8,   101.079200_r8,   100.113000_r8, &
                             74.076200_r8,    72.061400_r8,   149.118600_r8,   150.126000_r8,   136.228400_r8, &
                            185.234000_r8,   186.241400_r8,    92.136200_r8,   108.135600_r8,   173.140600_r8, &
                            174.148000_r8,   190.147400_r8,    98.098200_r8,    58.035600_r8,   121.047940_r8, &
                            119.074340_r8,   147.084740_r8,   162.117940_r8,   147.125940_r8,    64.064800_r8, &
                             62.132400_r8,    17.028940_r8,    18.036340_r8,   222.000000_r8,   207.200000_r8, &
                             27.025140_r8,    41.050940_r8,    26.036800_r8,    46.024600_r8,    63.031400_r8, &
                             98.078400_r8,    12.011000_r8,   144.132000_r8,   115.107340_r8,    12.011000_r8, &
                             12.011000_r8,    12.011000_r8,   135.064039_r8,    58.442468_r8,     1.007400_r8, &
                            115.107340_r8,    12.011000_r8,    58.442468_r8,     1.007400_r8,   135.064039_r8, &
                             58.442468_r8,   115.107340_r8,     1.007400_r8 /)

      crb_mass(:108) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                             12.011000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8, &
                             24.022000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8, &
                             24.022000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8,    36.033000_r8, &
                             36.033000_r8,    36.033000_r8,    36.033000_r8,    36.033000_r8,    36.033000_r8, &
                             36.033000_r8,    36.033000_r8,    36.033000_r8,    48.044000_r8,    48.044000_r8, &
                             48.044000_r8,    48.044000_r8,    48.044000_r8,    60.055000_r8,    60.055000_r8, &
                             60.055000_r8,    60.055000_r8,    60.055000_r8,    60.055000_r8,    48.044000_r8, &
                             48.044000_r8,    48.044000_r8,    48.044000_r8,    48.044000_r8,    60.055000_r8, &
                             36.033000_r8,    36.033000_r8,    60.055000_r8,    60.055000_r8,   120.110000_r8, &
                            120.110000_r8,   120.110000_r8,    84.077000_r8,    84.077000_r8,    84.077000_r8, &
                             84.077000_r8,    84.077000_r8,    60.055000_r8,    24.022000_r8,    24.022000_r8, &
                             36.033000_r8,    48.044000_r8,    60.055000_r8,    60.055000_r8,     0.000000_r8, &
                             24.022000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,    24.022000_r8,    24.022000_r8,    12.011000_r8,    12.011000_r8, &
                              0.000000_r8,    12.011000_r8,   144.132000_r8,     0.000000_r8,    12.011000_r8, &
                             12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8 /)

      fix_mass(:  4) = (/ 0.00000000_r8, 28.0134800_r8, 31.9988000_r8, 18.0142000_r8 /)

      clsmap(:  4,1) = (/   16,   4,  84,  85 /)
      clsmap(:104,4) = (/    1,   2,   3,  17,  11,  86,  87,  15,   5,   6, &
                             7,   8,   9,  10,  12,  13,  14,  18,  19,  20, &
                            21,  22,  28,  23,  24,  25,  26,  27,  29,  30, &
                            31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                            41,  42,  43,  44,  45,  49,  50,  51,  46,  47, &
                            48,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                            61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                            71,  72,  73,  74,  75,  76,  77,  78,  79,  80, &
                            81,  82,  83,  93,  88,  89,  90,  91,  92,  94, &
                            95,  96,  97,  98,  99, 100, 101, 102, 103, 104, &
                           105, 106, 107, 108 /)

      permute(:104,4) = (/   97,  90,  81, 104,  67,  32,  25,  78,  98,  99, &
                             96,  58,  46,  38, 100, 101,  37, 103,  48, 102, &
                             60,  35,  27,  45,  36,  66,  55,  75,  21,  74, &
                             43,  82,  95,  54,  80,  22,  76,  47,  71,  59, &
                             72,  85,  39,  23,  40,  24,  73,  69,  53,  70, &
                             41,  77,  93,  64,  91,  87,  92,  42,  94,  49, &
                             86,  88,  89,  33,  65,  79,  56,  28,  29,  61, &
                             50,  31,  62,  68,  57,  51,  63,  84,  83,  26, &
                             34,  20,   1,   2,  30,  44,  52,   3,   4,   5, &
                              6,   7,   8,   9,  10,  11,  12,  13,  14,  15, &
                             16,  17,  18,  19 /)

      diag_map(:104) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                            22,  25,  28,  31,  34,  38,  40,  43,  48,  51, &
                            57,  61,  65,  69,  74,  78,  82,  86,  92,  97, &
                           104, 109, 114, 119, 123, 130, 135, 140, 145, 148, &
                           154, 159, 165, 170, 175, 178, 185, 192, 196, 203, &
                           208, 216, 222, 230, 238, 249, 256, 261, 266, 277, &
                           285, 293, 303, 318, 327, 334, 345, 360, 368, 381, &
                           397, 406, 415, 425, 438, 448, 454, 465, 478, 492, &
                           506, 525, 545, 562, 587, 623, 657, 700, 723, 806, &
                           865, 881, 915, 928 /)

      extfrc_lst(:  9) = (/ 'CO              ','NO2             ','SO2             ','so4_a1          ','so4_a2          ', &
                            'pom_a1          ','bc_a1           ','num_a1          ','num_a2          ' /)

      frc_from_dataset(:  9) = (/ .true., .true., .true., .true., .true., &
                                  .true., .true., .true., .true. /)

      inv_lst(:  4) = (/ 'M               ', 'N2              ', 'O2              ', 'H2O             ' /)

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
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo2_b           ', 'jo3_a           ', 'jo3_b           ', 'jn2o            ', &
                                     'jno2            ', 'jn2o5_a         ', 'jn2o5_b         ', 'jhno3           ', &
                                     'jno3_a          ', 'jno3_b          ', 'jho2no2_a       ', 'jho2no2_b       ', &
                                     'jch3ooh         ', 'jch2o_a         ', 'jch2o_b         ', 'jh2o2           ', &
                                     'jch3cho         ', 'jpooh           ', 'jch3co3h        ', 'jpan            ', &
                                     'jmpan           ', 'jmacr_a         ', 'jmacr_b         ', 'jmvk            ', &
                                     'jc2h5ooh        ', 'jeooh           ', 'jc3h7ooh        ', 'jrooh           ', &
                                     'jacet           ', 'jmgly           ', 'jxooh           ', 'jonitr          ', &
                                     'jisopooh        ', 'jhyac           ', 'jglyald         ', 'jmek            ', &
                                     'jbigald         ', 'jglyoxal        ', 'jalkooh         ', 'jmekooh         ', &
                                     'jtolooh         ', 'jterpooh        ', 'usr_O_O2        ', 'O_O3            ', &
                                     'O1D_N2          ', 'O1D_O2b         ', 'O1D_O3          ', 'ox_l1           ', &
                                     'O1D_H2          ', 'O1D_N2Oa        ', 'O1D_N2Ob        ', 'O1D_CH4a        ', &
                                     'O1D_CH4b        ', 'O1D_CH4c        ', 'O1D_HCN         ', 'H_O2            ', &
                                     'H_O3            ', 'H_HO2a          ', 'H_HO2b          ', 'H_HO2c          ', &
                                     'OH_O            ', 'ox_l2           ', 'OH_HO2          ', 'OH_OH           ', &
                                     'OH_OH_M         ', 'OH_H2           ', 'OH_H2O2         ', 'H2_O            ', &
                                     'HO2_O           ', 'ox_l3           ', 'usr_HO2_HO2     ', 'H2O2_O          ', &
                                     'HCN_OH          ', 'CH3CN_OH        ', 'NO_O_M          ', 'ox_p1           ', &
                                     'NO_O3           ', 'NO2_O           ', 'NO2_O_M         ', 'NO2_O3          ', &
                                     'tag_NO2_NO3     ', 'usr_N2O5_M      ', 'tag_NO2_OH      ', 'usr_HNO3_OH     ', &
                                     'NO3_NO          ', 'NO3_O           ', 'NO3_OH          ', 'NO3_HO2         ', &
                                     'tag_NO2_HO2     ', 'HO2NO2_OH       ', 'usr_HO2NO2_M    ', 'CH4_OH          ', &
                                     'usr_CO_OH_b     ', 'CO_OH_M         ', 'CH2O_NO3        ', 'CH2O_OH         ', &
                                     'CH2O_O          ', 'CH2O_HO2        ', 'ox_p2           ', 'CH3O2_HO2       ', &
                                     'CH3O2_CH3O2a    ', 'CH3O2_CH3O2b    ', 'CH3OH_OH        ', 'CH3OOH_OH       ', &
                                     'HCOOH_OH        ', 'HOCH2OO_M       ', 'HOCH2OO_NO      ', 'HOCH2OO_HO2     ', &
                                     'C2H2_OH_M       ', 'C2H6_OH         ', 'tag_C2H4_OH     ', 'ox_p16          ', &
                                     'EO2_HO2         ', 'EO_O2           ', 'EO_M            ', 'ox_l6           ', &
                                     'CH3COOH_OH      ', 'ox_p5           ', 'C2H5O2_HO2      ', 'C2H5O2_CH3O2    ', &
                                     'C2H5O2_C2H5O2   ', 'C2H5OOH_OH      ', 'CH3CHO_OH       ', 'CH3CHO_NO3      ', &
                                     'ox_p4           ', 'tag_CH3CO3_NO2  ', 'CH3CO3_HO2      ', 'CH3CO3_CH3O2    ', &
                                     'CH3CO3_CH3CO3   ', 'CH3COOOH_OH     ', 'GLYALD_OH       ', 'GLYOXAL_OH      ', &
                                     'C2H5OH_OH       ', 'usr_PAN_M       ', 'PAN_OH          ', 'tag_C3H6_OH     ', &
                                     'ox_l4           ', 'C3H6_NO3        ', 'ox_p9           ', 'C3H7O2_HO2      ', &
                                     'CH3H7O2_CH3O2   ', 'CH3H7OOH_OH     ', 'C3H8_OH         ', 'ox_p3           ', &
                                     'PO2_HO2         ', 'POOH_OH         ', 'usr_CH3COCH3_OH ', 'ox_p10          ', &
                                     'RO2_HO2         ', 'RO2_CH3O2       ', 'ROOH_OH         ', 'HYAC_OH         ', &
                                     'CH3COCHO_OH     ', 'CH3COCHO_NO3    ', 'ONIT_OH         ', 'BIGENE_OH       ', &
                                     'ox_p15          ', 'MVK_OH          ', 'ox_l7           ', 'MEK_OH          ', &
                                     'ox_p17          ', 'MEKO2_HO2       ', 'MEKOOH_OH       ', 'MACR_OH         ', &
                                     'ox_l8           ', 'ox_p7           ', 'MACRO2_NOb      ', 'MACRO2_NO3      ', &
                                     'MACRO2_HO2      ', 'MACRO2_CH3O2    ', 'MACRO2_CH3CO3   ', 'MACROOH_OH      ', &
                                     'ox_p8           ', 'MCO3_NO3        ', 'MCO3_HO2        ', 'MCO3_CH3O2      ', &
                                     'MCO3_CH3CO3     ', 'MCO3_MCO3       ', 'usr_MCO3_NO2    ', 'usr_MPAN_M      ', &
                                     'MPAN_OH_M       ', 'ISOP_OH         ', 'ox_l5           ', 'ISOP_NO3        ', &
                                     'ox_p6           ', 'ISOPO2_NO3      ', 'ISOPO2_HO2      ', 'ISOPOOH_OH      ', &
                                     'ISOPO2_CH3O2    ', 'ISOPO2_CH3CO3   ', 'ISOPNO3_NO      ', 'ISOPNO3_NO3     ', &
                                     'ISOPNO3_HO2     ', 'BIGALK_OH       ', 'ONITR_OH        ', 'ONITR_NO3       ', &
                                     'HYDRALD_OH      ', 'ox_p14          ', 'ALKO2_HO2       ', 'ALKOOH_OH       ', &
                                     'ox_p11          ', 'XO2_NO3         ', 'XO2_HO2         ', 'XO2_CH3O2       ', &
                                     'XO2_CH3CO3      ', 'XOOH_OHa        ', 'usr_XOOH_OH     ', 'TOLUENE_OH      ', &
                                     'ox_p12          ', 'TOLO2_HO2       ', 'TOLO2_OH        ', 'CRESOL_OH       ', &
                                     'XOH_NO2         ', 'C10H16_OH       ', 'soa1            ', 'C10H16_NO3      ', &
                                     'ox_p13          ', 'TERPO2_HO2      ', 'TERPOOH_OH      ', 'usr_N2O5_aer    ', &
                                     'usr_NO3_aer     ', 'usr_NO2_aer     ', 'usr_SO2_OH      ', 'DMS_OHb         ', &
                                     'usr_DMS_OH      ', 'DMS_NO3         ', 'NH3_OH          ', 'usr_HO2_aer     ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                                       11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                                       21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                                       31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                                       41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                                       51,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                                       61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                                       71,  72,  73,  74,  75,  76,  77,  78,  79,  80, &
                                       81,  82,  83,  84,  85,  86,  87,  88,  89,  90, &
                                       91,  92,  93,  94,  95,  96,  97,  98,  99, 100, &
                                      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, &
                                      111, 112, 113, 114, 115, 116, 117, 118, 119, 120, &
                                      121, 122, 123, 124, 125, 126, 127, 128, 129, 130, &
                                      131, 132, 133, 134, 135, 136, 137, 138, 139, 140, &
                                      141, 142, 143, 144, 145, 146, 147, 148, 149, 150, &
                                      151, 152, 153, 154, 155, 156, 157, 158, 159, 160, &
                                      161, 162, 163, 164, 165, 166, 167, 168, 169, 170, &
                                      171, 172, 173, 174, 175, 176, 177, 178, 179, 180, &
                                      181, 182, 183, 184, 185, 186, 187, 188, 189, 190, &
                                      191, 192, 193, 194, 195, 196, 197, 198, 199, 200, &
                                      201, 202, 203, 204, 205, 206, 207, 208, 209, 210, &
                                      211, 212, 213, 214, 215, 216, 217, 218, 219, 221, &
                                      222, 223, 224, 225, 226, 227, 228, 229 /)
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
      pht_alias_lst(:,1) = (/ 'userdefined     ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ' /)
      pht_alias_lst(:,2) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', 'jch3ooh         ', 'jh2o2           ', '                ', &
                              'jpan            ', '                ', '                ', '                ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', &
                              '                ', '                ', 'jch3ooh         ', 'jch3cho         ', &
                              'jch3ooh         ', '                ', '                ', 'jacet           ', &
                              'jno2            ', 'jmgly           ', 'jch3ooh         ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jch3ooh         ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, .28_r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, .2_r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
