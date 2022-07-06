
      module chem_mods
!--------------------------------------------------------------
!     	... Basic chemistry parameters and arrays
!--------------------------------------------------------------

      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

      save

      integer, parameter :: phtcnt    = PHTCNT, &    ! number of photolysis reactions
                            rxntot    = RXNCNT, &    ! number of total reactions
                            gascnt    = GASCNT, &    ! number of gas phase reactions
                            nabscol   = NCOL, &      ! number of absorbing column densities
                            gas_pcnst = PCNST, &     ! number of "gas phase" species
                            nfs       = NFS, &       ! number of "fixed" species
                            relcnt    = RELCNT, &    ! number of relationship species
                            grpcnt    = GRPCNT, &    ! number of group members
                            nzcnt     = IMP_NZCNT, & ! number of non-zero matrix entries
                            extcnt    = EXTCNT, &    ! number of species with external forcing
                            clscnt1   = CLSCNT1, &   ! number of species in explicit class
                            clscnt2   = CLSCNT2, &   ! number of species in hov class
                            clscnt3   = CLSCNT3, &   ! number of species in ebi class
                            clscnt4   = CLSCNT4, &   ! number of species in implicit class
                            clscnt5   = CLSCNT5, &   ! number of species in rodas class
                            indexm    = INDEXM, &    ! index of total atm density in invariant array
                            indexh2o  = INDEXH2O, &  ! index of water vapor density
                            clsze     = CLSZE, &     ! loop length for implicit chemistry
                            rxt_tag_cnt = RXTTAGCNT, &
                            nslvd     = NSLVD

      integer   :: clscnt(5)            = 0
      integer   :: cls_rxt_cnt(4,5)     = 0
      integer   :: clsmap(gas_pcnst,5)  = 0
      integer   :: permute(gas_pcnst,5) = 0
# if CLSCNT4 != 0
      integer   :: diag_map(clscnt4)    = 0
# elif CLSCNT5 != 0
      integer   :: diag_map(clscnt5)    = 0
# endif
      real(r8)  :: adv_mass(gas_pcnst)  = 0._r8
      real(r8)  :: crb_mass(gas_pcnst)  = 0._r8
      real(r8)  :: fix_mass(max(1,nfs))
# if GRPCNT != 0
      real(r8)  :: nadv_mass(grpcnt)    = 0._r8
# endif

      integer, allocatable :: rxt_tag_map(:)
      real(r8), allocatable :: pht_alias_mult(:,:)
      character(len=16), allocatable :: rxt_tag_lst(:)
      character(len=16), allocatable :: pht_alias_lst(:,:)
      character(len=16)              :: inv_lst(max(1,nfs))
      character(len=16)              :: extfrc_lst(max(1,extcnt))
      logical                        :: frc_from_dataset(max(1,extcnt))
      character(len=16)              :: slvd_lst(max(1,nslvd))

      end module chem_mods
