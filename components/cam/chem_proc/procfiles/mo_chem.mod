
      module chem_mods
!--------------------------------------------------------------
!     	... basic chemistry array parameters
!--------------------------------------------------------------

      use mo_grid, only : pcnstm1

      implicit none

      save

      integer, parameter :: hetcnt     = HETCNT, &    ! number of heterogeneous processes
                            phtcnt     = PHTCNT, &    ! number of photo processes
                            rxntot     = RXNCNT, &    ! number of total reactions
                            gascnt     = GASCNT, &    ! number of gas phase reactions
                            nfs        = NFS, &       ! number of "fixed" species
                            relcnt     = RELCNT, &    ! number of relationship species
                            grpcnt     = GRPCNT, &    ! number of group members
                            imp_nzcnt  = IMP_NZCNT, &     ! number of non-zero implicit matrix entries
                            rod_nzcnt  = ROD_NZCNT, &     ! number of non-zero rodas matrix entries
                            extcnt     = EXTCNT, &    ! number of species with external forcing
                            clscnt1    = CLSCNT1, &  ! number of species in explicit class
                            clscnt2    = CLSCNT2, &  ! number of species in hov class
                            clscnt3    = CLSCNT3, &  ! number of species in ebi class
                            clscnt4    = CLSCNT4, &  ! number of species in implicit class
                            clscnt5    = CLSCNT5, &  ! number of species in rodas class
                            indexm     = INDEXM, &    ! index of total atm density in invariant array
                            ncol_abs   = NCOL, &    ! number of column densities
                            indexh2o   = INDEXH2O, &    ! index of water vapor density
                            clsze      = CLSZE       ! loop length for implicit chemistry

      integer ::            ngrp       = 0
      integer ::            drydep_cnt = 0
      integer ::            srfems_cnt = 0
      integer ::            rxt_alias_cnt = 0
      integer ::            fbc_cnt(2) = 0
      integer, allocatable :: grp_mem_cnt(:)
      integer, allocatable :: rxt_alias_map(:)
      real      :: adv_mass(pcnstm1)
      real      :: nadv_mass(grpcnt)
      character(len=16), allocatable :: rxt_alias_lst(:)
      character(len=8), allocatable  :: drydep_lst(:)
      character(len=8), allocatable  :: srfems_lst(:)
      character(len=8), allocatable  :: grp_lst(:)
      character(len=8), allocatable  :: flbc_lst(:)
      character(len=8)               :: het_lst(max(1,hetcnt))
      character(len=8)               :: extfrc_lst(max(1,extcnt))
      character(len=8)               :: inv_lst(max(1,nfs))

      type solver_class
	 integer :: clscnt
	 integer :: lin_rxt_cnt
	 integer :: nln_rxt_cnt
	 integer :: indprd_cnt
	 integer :: iter_max
         integer :: cls_rxt_cnt(4)
         integer, pointer :: permute(:)
         integer, pointer :: diag_map(:)
         integer, pointer :: clsmap(:)
      end type solver_class

      type(solver_class) :: explicit, implicit, rodas

      contains

      subroutine chem_mods_inti
!--------------------------------------------------------------
!     	... intialize the class derived type
!--------------------------------------------------------------

      implicit none

      integer :: astat

      explicit%clscnt       = CLSCNT1
      explicit%indprd_cnt   = CLSINDPRD1

      implicit%clscnt       = CLSCNT4
      implicit%lin_rxt_cnt  = IMP_LINCNT
      implicit%nln_rxt_cnt  = IMP_NLNCNT
      implicit%indprd_cnt   = CLSINDPRD4
      implicit%iter_max     = IMPITERMAX

      rodas%clscnt          = CLSCNT5
      rodas%lin_rxt_cnt     = ROD_LINCNT
      rodas%nln_rxt_cnt     = ROD_NLNCNT
      rodas%indprd_cnt      = CLSINDPRD5

      if( explicit%clscnt > 0 ) then
	 allocate( explicit%clsmap(explicit%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_inti: failed to allocate explicit%clsmap ; error = ',astat
	    call endrun
	 end if
         explicit%clsmap(:)  = 0
      end if
      if( implicit%clscnt > 0 ) then
	 allocate( implicit%permute(implicit%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_inti: failed to allocate implicit%permute ; error = ',astat
	    call endrun
	 end if
         implicit%permute(:)  = 0
	 allocate( implicit%diag_map(implicit%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_inti: failed to allocate implicit%diag_map ; error = ',astat
	    call endrun
	 end if
         implicit%diag_map(:)  = 0
	 allocate( implicit%clsmap(implicit%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_inti: failed to allocate implicit%clsmap ; error = ',astat
	    call endrun
	 end if
         implicit%clsmap(:)  = 0
      end if
      if( rodas%clscnt > 0 ) then
	 allocate( rodas%permute(rodas%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_inti: failed to allocate rodas%permute ; error = ',astat
	    call endrun
	 end if
         rodas%permute(:)  = 0
	 allocate( rodas%diag_map(rodas%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_inti: failed to allocate rodas%diag_map ; error = ',astat
	    call endrun
	 end if
         rodas%diag_map(:)  = 0
	 allocate( rodas%clsmap(rodas%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_inti: failed to allocate rodas%clsmap ; error = ',astat
	    call endrun
	 end if
         rodas%clsmap(:)  = 0
      end if

      end subroutine chem_mods_inti

      end module chem_mods
