













      module chem_mods
!--------------------------------------------------------------
!     	... Basic chemistry parameters and arrays
!--------------------------------------------------------------

      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

      save

      integer, parameter :: phtcnt    = 1, &    ! number of photolysis reactions
                            rxntot    = 7, &    ! number of total reactions
                            gascnt    = 6, &    ! number of gas phase reactions
                            nabscol   = 2, &      ! number of absorbing column densities
                            gas_pcnst = 56, &     ! number of "gas phase" species
                            nfs       = 8, &       ! number of "fixed" species
                            relcnt    = 0, &    ! number of relationship species
                            grpcnt    = 0, &    ! number of group members
                            nzcnt     = 57, & ! number of non-zero matrix entries
                            extcnt    = 9, &    ! number of species with external forcing
                            clscnt1   = 1, &   ! number of species in explicit class
                            clscnt2   = 0, &   ! number of species in hov class
                            clscnt3   = 0, &   ! number of species in ebi class
                            clscnt4   = 55, &   ! number of species in implicit class
                            clscnt5   = 0, &   ! number of species in rodas class
                            indexm    = 1, &    ! index of total atm density in invariant array
                            indexh2o  = 4, &  ! index of water vapor density
                            clsze     = 1, &     ! loop length for implicit chemistry
                            rxt_tag_cnt = 4, &
                            nslvd     = 0

      integer   :: clscnt(5)            = 0
      integer   :: cls_rxt_cnt(4,5)     = 0
      integer   :: clsmap(gas_pcnst,5)  = 0
      integer   :: permute(gas_pcnst,5) = 0
      integer   :: diag_map(clscnt4)    = 0
      real(r8)  :: adv_mass(gas_pcnst)  = 0._r8
      real(r8)  :: crb_mass(gas_pcnst)  = 0._r8
      real(r8)  :: fix_mass(max(1,nfs))

      integer, allocatable :: rxt_tag_map(:)
      real(r8), allocatable :: pht_alias_mult(:,:)
      character(len=16), allocatable :: rxt_tag_lst(:)
      character(len=16), allocatable :: pht_alias_lst(:,:)
      character(len=16)              :: inv_lst(max(1,nfs))
      character(len=16)              :: extfrc_lst(max(1,extcnt))
      logical                        :: frc_from_dataset(max(1,extcnt))
      character(len=16)              :: slvd_lst(max(1,nslvd))

      end module chem_mods
  
      module m_spc_id                                                           
  
      implicit none                                                             
  
      integer, parameter :: id_O3 =   1
      integer, parameter :: id_H2O2 =   2
      integer, parameter :: id_H2SO4 =   3
      integer, parameter :: id_SO2 =   4
      integer, parameter :: id_DMS =   5
      integer, parameter :: id_SOAG =   6
      integer, parameter :: id_so4_a1 =   7
      integer, parameter :: id_pom_a1 =   8
      integer, parameter :: id_soa_a1 =   9
      integer, parameter :: id_bc_a1 =  10
      integer, parameter :: id_dst_a1 =  11
      integer, parameter :: id_ncl_a1 =  12
      integer, parameter :: id_mom_a1 =  13
      integer, parameter :: id_num_a1 =  14
      integer, parameter :: id_so4_a2 =  15
      integer, parameter :: id_soa_a2 =  16
      integer, parameter :: id_ncl_a2 =  17
      integer, parameter :: id_mom_a2 =  18
      integer, parameter :: id_num_a2 =  19
      integer, parameter :: id_dst_a3 =  20
      integer, parameter :: id_ncl_a3 =  21
      integer, parameter :: id_so4_a3 =  22
      integer, parameter :: id_bc_a3 =  23
      integer, parameter :: id_pom_a3 =  24
      integer, parameter :: id_soa_a3 =  25
      integer, parameter :: id_mom_a3 =  26
      integer, parameter :: id_num_a3 =  27
      integer, parameter :: id_pom_a4 =  28
      integer, parameter :: id_bc_a4 =  29
      integer, parameter :: id_mom_a4 =  30
      integer, parameter :: id_num_a4 =  31
      integer, parameter :: id_so4_c1 =  32
      integer, parameter :: id_pom_c1 =  33
      integer, parameter :: id_soa_c1 =  34
      integer, parameter :: id_bc_c1 =  35
      integer, parameter :: id_dst_c1 =  36
      integer, parameter :: id_ncl_c1 =  37
      integer, parameter :: id_mom_c1 =  38
      integer, parameter :: id_num_c1 =  39
      integer, parameter :: id_so4_c2 =  40
      integer, parameter :: id_soa_c2 =  41
      integer, parameter :: id_ncl_c2 =  42
      integer, parameter :: id_mom_c2 =  43
      integer, parameter :: id_num_c2 =  44
      integer, parameter :: id_dst_c3 =  45
      integer, parameter :: id_ncl_c3 =  46
      integer, parameter :: id_so4_c3 =  47
      integer, parameter :: id_bc_c3 =  48
      integer, parameter :: id_pom_c3 =  49
      integer, parameter :: id_soa_c3 =  50
      integer, parameter :: id_mom_c3 =  51
      integer, parameter :: id_num_c3 =  52
      integer, parameter :: id_pom_c4 =  53
      integer, parameter :: id_bc_c4 =  54
      integer, parameter :: id_mom_c4 =  55
      integer, parameter :: id_num_c4 =  56
  
  
      end module m_spc_id                                                       
                                                                                
      module m_rxt_id                                                           
                                                                                
      implicit none                                                             
                                                                                
      integer, parameter :: rid_jh2o2 =    1                                    
      integer, parameter :: rid_usr_HO2_HO2 =    2                              
      integer, parameter :: rid_usr_SO2_OH =    4                               
      integer, parameter :: rid_usr_DMS_OH =    6                               
                                                                                
      integer, parameter :: rid_r0003 =    3                                    
      integer, parameter :: rid_r0005 =    5                                    
      integer, parameter :: rid_r0007 =    7                                    
                                                                                
      end module m_rxt_id                                                       
