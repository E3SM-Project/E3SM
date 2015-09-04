
      module io

      implicit none

      integer :: lin      ! input unit number
      integer :: lout     ! output unit number

      character(len=320) :: buff        ! primary line input buffer
      character(len=320) :: buffh       ! upcase xform of buff
      character(len=320) :: procout_path   = "../output/"
      character(len=320) :: procfiles_path = "../procfiles/cam/"
      character(len=320) :: output_path    = "../output/"
      character(len=320) :: input_path
      character(len=320) :: temp_path
      character(len=320) :: src_dir      = "../bkend/"
      character(len=320) :: sim_dat_path = "../output/"
      character(len=320) :: sim_dat_filespec
      character(len=320) :: sim_dat_filename

      end module io

!-----------------------------------------------------------
!	... Table of the elements; symbol and amu
!-----------------------------------------------------------
      module elements

      integer, private,parameter   :: dp = selected_real_kind( 12 )

      type element
         character(len=2) :: sym
         real(dp)         :: wght
      end type element

      integer, private             :: tab_max = 100
      integer, private             :: id_cnt = 1
      character(len=39), private   :: id
      type( element ), private     :: e_table(100)

      contains

      subroutine iniele()
!-----------------------------------------------------------
!	... Initialize the element mass table and mass computation
!-----------------------------------------------------------

      implicit none

!-----------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------
      integer :: i

      e_table(:)%sym = '  '
      e_table(1) =  ELEMENT( 'H ',1.0074_dp )
      e_table(2) =  ELEMENT( 'He',4.0020602_dp )
      e_table(3) =  ELEMENT( 'Li',6.941_dp )
      e_table(4) =  ELEMENT( 'Be',9.012182_dp )
      e_table(5) =  ELEMENT( 'B ',10.811_dp )
      e_table(6) =  ELEMENT( 'C ',12.011_dp )
      e_table(7) =  ELEMENT( 'N ',14.00674_dp )
      e_table(8) =  ELEMENT( 'O ',15.9994_dp )
      e_table(9) =  ELEMENT( 'F ',18.9984032_dp )
      e_table(10) = ELEMENT( 'Ne',20.1797_dp )
      e_table(11) = ELEMENT( 'Na',22.989768_dp )
      e_table(12) = ELEMENT( 'Mg',24.305_dp )
      e_table(13) = ELEMENT( 'Al',26.981539_dp )
      e_table(14) = ELEMENT( 'Si',28.0855_dp )
      e_table(15) = ELEMENT( 'P ',30.97362_dp )
      e_table(16) = ELEMENT( 'S ',32.066_dp )
      e_table(17) = ELEMENT( 'Cl',35.4527_dp )
      e_table(18) = ELEMENT( 'Ar',39.948_dp )
      e_table(19) = ELEMENT( 'K ',39.0983_dp )
      e_table(20) = ELEMENT( 'Ca',40.078_dp )
      e_table(21) = ELEMENT( 'Sc',44.95591_dp )
      e_table(22) = ELEMENT( 'Ti',47.867_dp )
      e_table(23) = ELEMENT( 'V ',50.9415_dp )
      e_table(24) = ELEMENT( 'Cr',51.9961_dp )
      e_table(25) = ELEMENT( 'Mn',54.93085_dp )
      e_table(26) = ELEMENT( 'Fe',55.845_dp )
      e_table(27) = ELEMENT( 'Co',58.9332_dp )
      e_table(28) = ELEMENT( 'Ni',58.6934_dp )
      e_table(29) = ELEMENT( 'Cu',63.546_dp )
      e_table(30) = ELEMENT( 'Zn',65.39_dp )
      e_table(31) = ELEMENT( 'Ga',69.723_dp )
      e_table(32) = ELEMENT( 'Ge',72.61_dp )
      e_table(33) = ELEMENT( 'As',74.92159_dp )
      e_table(34) = ELEMENT( 'Se',78.96_dp )
      e_table(35) = ELEMENT( 'Br',79.904_dp )
      e_table(36) = ELEMENT( 'Kr',83.8_dp )
      e_table(37) = ELEMENT( 'Rb',85.4678_dp )
      e_table(38) = ELEMENT( 'Sr',87.62_dp )
      e_table(39) = ELEMENT( 'Y ',88.90585_dp )
      e_table(40) = ELEMENT( 'Zr',91.224_dp )
      e_table(41) = ELEMENT( 'Nb',92.90638_dp )
      e_table(42) = ELEMENT( 'Mo',95.94_dp )
      e_table(43) = ELEMENT( 'Tc',98._dp )
      e_table(44) = ELEMENT( 'Ru',101.07_dp )
      e_table(45) = ELEMENT( 'Rh',102.9055_dp )
      e_table(46) = ELEMENT( 'Pd',106.42_dp )
      e_table(47) = ELEMENT( 'Ag',107.8682_dp )
      e_table(48) = ELEMENT( 'Cd',112.411_dp )
      e_table(49) = ELEMENT( 'In',114.818_dp )
      e_table(50) = ELEMENT( 'Sn',118.71_dp )
      e_table(51) = ELEMENT( 'Sb',121.76_dp )
      e_table(52) = ELEMENT( 'Te',127.6_dp )
      e_table(53) = ELEMENT( 'I ',126.90447_dp )
      e_table(54) = ELEMENT( 'Xe',131.29_dp )
      e_table(55) = ELEMENT( 'Cs',132.90543_dp )
      e_table(56) = ELEMENT( 'Ba',137.327_dp )
      e_table(57) = ELEMENT( 'La',138.9055_dp )
      e_table(58) = ELEMENT( 'Hf',178.49_dp )
      e_table(59) = ELEMENT( 'Ta',180.9479_dp )
      e_table(60) = ELEMENT( 'W ',183.84_dp )
      e_table(61) = ELEMENT( 'Re',186.207_dp )
      e_table(62) = ELEMENT( 'Os',190.23_dp )
      e_table(63) = ELEMENT( 'Ir',192.217_dp )
      e_table(64) = ELEMENT( 'Pt',195.08_dp )
      e_table(65) = ELEMENT( 'Au',196.96654_dp )
      e_table(66) = ELEMENT( 'Hg',200.59_dp )
      e_table(67) = ELEMENT( 'Tl',204.3833_dp )
      e_table(68) = ELEMENT( 'Pb',207.2_dp )
      e_table(69) = ELEMENT( 'Bi',208.98037_dp )
      e_table(70) = ELEMENT( 'Po',209._dp )
      e_table(71) = ELEMENT( 'At',210._dp )
      e_table(72) = ELEMENT( 'Rn',222._dp )
      e_table(73) = ELEMENT( 'Fr',223._dp )
      e_table(74) = ELEMENT( 'Ra',226.025_dp )
      e_table(75) = ELEMENT( 'Ac',227.028_dp )
      e_table(75) = ELEMENT( 'e',.000548567_dp )

      do i = 1,100
	 if( e_table(i)%sym == '  ' ) then
	    exit
	 end if
      end do
      tab_max = i - 1

      id(:1) = e_table(1)%sym(:1)
      do i = 2,tab_max
	 if( scan( e_table(i)%sym(:1), id(:id_cnt) ) == 0 ) then
	    id_cnt = id_cnt + 1
	    id(id_cnt:id_cnt) = e_table(i)%sym(:1)
	 end if
      end do

      end subroutine iniele

      real(dp) function com_mass( compound, carbmass )
!-----------------------------------------------------------
!	... Compute the mass of input compound
!-----------------------------------------------------------

      implicit none

!-----------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------
      character(len=*), intent(in) :: compound
      logical, optional,intent(in) :: carbmass

!-----------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------
      integer  :: beg, end, pos, nump, index
      integer  :: ios, el_cnt
      real(dp) :: sum
      logical :: carbon_only
      
      carbon_only = .false.
      if (present(carbmass)) then
         carbon_only = carbmass
      endif

      end = len_trim( compound )
      sum = 0._dp
      do
	 pos = scan( compound(:end), id(:id_cnt), back = .true. )
	 if( pos == 0 ) then
	    exit
	 end if
	 nump = scan( compound(pos+1:end), '0123456789' )
	 if( nump /= 0 ) then
	    nump = pos + nump
	    read(compound(nump:end),*,iostat=ios) el_cnt
	    if( ios /= 0 .or. el_cnt == 0 ) then
	       exit
	    end if
	    end = nump - 1
	 else
	    el_cnt = 1
	 end if
	 do index = 1,tab_max
	    if( e_table(index)%sym == compound(pos:end) ) then
	       exit
	    end if
	 end do
!        if( index > 39 ) then
!           COM_MASS = 0.
!           exit
!        end if
         if (carbon_only) then
           if( trim(e_table(index)%sym) == 'C') then
	     sum = sum + e_table(index)%wght * real( el_cnt,dp )
           endif
         else
	   sum = sum + e_table(index)%wght * real( el_cnt,dp )
         endif
	 end = pos - 1
	 if( end <= 0 ) then
	    exit
	 end if
      end do
      com_mass = sum

      end function com_mass

      end module elements

      module SP_MODS
	 integer :: n  = 0        ! order of matrix
	 integer :: nz = 0        ! # of non=zero elements
	 integer :: sp = 0        ! stack pointer
	 integer :: nb = 0        ! search counter
	 integer :: pp = 0        ! perm vector position
	 integer :: blkcnt = 0    ! strongly connected blk count
	 integer, allocatable :: number(:)
	 integer, allocatable :: lowlink(:)
	 integer, allocatable :: vstack(:)
	 integer, allocatable :: perm(:)
	 integer, allocatable :: rp(:)
	 integer, allocatable :: ci(:)
	 integer, allocatable :: stcoblk(:)
	 integer, allocatable :: blkmemcnt(:)
	 logical, allocatable :: matrix(:,:)

         type SPARSITY
            integer, pointer :: diag_map(:)                 ! map of jacobian diagonals
            integer, pointer :: mat_sp_map(:,:)             ! matrix sparsity "map"
            logical, pointer, dimension(:,:) :: mat_sp_pat, lu_sp_pat
         end type SPARSITY
      end module SP_MODS

      module VAR_MOD
!-----------------------------------------------------------------------
!        ... Mozart reaction variables
!-----------------------------------------------------------------------

      implicit none
 
      integer, parameter   :: dp = selected_real_kind( 12 )

      integer, parameter :: var_lim = 1000
      integer :: hst_file_lim
      integer :: hst_map_lim
      integer, pointer :: nq, relcnt, nfs, ngrp, &
			  ncol, new_nq, grp_mem_cnt 
      integer, target :: spccnt(7) = 0
      integer, allocatable :: &
                   grpflg(:), &
                   mem2grp_map(:), &
                   newind(:), &
                   grpmap(:,:), &
                   relmap(:,:), &
                   grpcnt(:), &
                   colmap(:), &
                   rel_flg(:), &
		   grp_rat_ind(:)

!-----------------------------------------------------------------------
!        ... The solution class variables
!-----------------------------------------------------------------------
      integer, allocatable :: &
                   clsmap(:,:,:)
      integer  ::  cls_ind_prdcnt
      integer  ::  clscnt(5) = 0                    ! count of solution species in each numerical "class"

      integer  ::  ptplen = 0                       ! total hist tape fields
      integer, allocatable :: &
                   histout_cnt(:,:,:), &
                   histout_map(:,:,:,:)
      integer, dimension(5,2) ::  &
		   class_prod_cnt = 0, &
		   class_loss_cnt = 0
      integer  ::  indexm = 0, &                    ! index for fixed species denoting total atm density
                   indexh2o = 0                     ! index for fixed species denoting water vapor density
      
      integer :: extcnt(5) = 0
      integer, allocatable :: &
                   srf_flx_map(:)                   ! surface flux flag
      integer  ::  srf_flx_cnt = 0                  ! count of soln species with surface emissions
      integer, allocatable :: &
                   dvel_map(:)                      ! deposition flux flag
      integer  ::  dvel_cnt = 0                     ! count of soln species with deposition flux

      real, allocatable :: &
                   colub(:), &                      ! upper boundary column integral
                   grpcof(:,:)                      ! multiplier for group members

      real(dp), allocatable :: &
                   mass(:), &                       ! molecular mass of the mechanism compound
                   c_mass(:), &                     ! carbon mass of the mechanism compound
                   temp_mass(:)                     ! original species masses and temp space

      character(len=64), allocatable :: aliases(:)
      character(len=16), target, allocatable  ::  spcsym(:,:)
      character(len=16), pointer ::  pcesym(:), &
                                     solsym(:), &
                                     fixsym(:), &
                                     grpsym(:), &
                                     colsym(:), &
                                     new_solsym(:), &
                                     grp_mem_sym(:), &
                                     slvdsym(:)

      integer :: nslvd

      character(len=16), allocatable :: &
                   user_hst_names(:,:)

      integer, allocatable :: permute(:,:)         ! permutation vector
      integer, allocatable :: permutation(:), permute_orig(:,:)

      contains

      subroutine VAR_INI()
!-----------------------------------------------------------------------
!        ... Allocate and initialize reaction variables
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer :: astat

      hst_file_lim = 10
      hst_map_lim  = 1000
      allocate( grpflg(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate grpflg'
	 stop
      end if
      grpflg(:) = 0
      allocate( mem2grp_map(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate mem2grp'
	 stop
      end if
      mem2grp_map(:) = 0
      allocate( newind(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate newind'
	 stop
      end if
      newind(:) = 0
      allocate( relmap(var_lim,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate relmap'
	 stop
      end if
      relmap(:,:) = 0
      allocate( grpcnt(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate grpcnt'
	 stop
      end if
      grpcnt(:) = 0
      allocate( grpmap(var_lim,var_lim/2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate grpmap'
	 stop
      end if
      grpmap(:,:) = 0
      allocate( colmap(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate colmap'
	 stop
      end if
      colmap(:) = 0
      allocate( rel_flg(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate rel_flg'
	 stop
      end if
      rel_flg(:) = 0
      allocate( grp_rat_ind(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate grp_rat_ind'
	 stop
      end if
      grp_rat_ind(:) = 0
      allocate( clsmap(var_lim,5,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate clsmap'
	 stop
      end if
      clsmap(:,:,:) = 0
      allocate( histout_cnt(20,2,hst_file_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate histout_cnt'
	 stop
      end if
      histout_cnt(:,:,:) = 0
      allocate( histout_map(hst_map_lim,20,2,hst_file_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate histout_map'
	 stop
      end if
      histout_map(:,:,:,:) = 0
      allocate( srf_flx_map(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate srf_flx_map'
	 stop
      end if
      srf_flx_map(:) = 0
      allocate( dvel_map(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate dvel_map'
	 stop
      end if
      dvel_map(:) = 0
      allocate( colub(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate colub'
	 stop
      end if
      colub(:) = 0.
      allocate( grpcof(var_lim,var_lim/2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate grpcof'
	 stop
      end if
      grpcof(:,:) = 1.
      allocate( mass(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate mass'
	 stop
      end if
      mass(:) = 0.
      allocate( c_mass(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate c_mass'
	 stop
      end if
      c_mass(:) = 0.
      allocate( temp_mass(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate temp_mass'
	 stop
      end if
      temp_mass(:) = 0.
      allocate( aliases(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate aliases'
	 stop
      end if
      aliases(:) = ' '
      allocate( spcsym(var_lim,8),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate spcsym'
	 stop
      end if
      spcsym(:,:) = ' '
      allocate( user_hst_names(var_lim,4),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate usr_hst_names'
	 stop
      end if
      user_hst_names(:,:) = ' '
      allocate( permute(var_lim,5),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate permute'
	 stop
      end if
      permute(:,:) = 0
      allocate( permute_orig(var_lim,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate permute_orig'
	 stop
      end if
      permute_orig(:,:) = 0
      allocate( permutation(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'VAR_INI: Failed to allocate permutation'
	 stop
      end if
      permutation(:) = 0

      nq     => spccnt(1)
      relcnt => spccnt(2)
      nfs    => spccnt(3)
      ngrp   => spccnt(4)
      ncol   => spccnt(5)
      new_nq => spccnt(6)
      grp_mem_cnt => spccnt(7)
      solsym => spcsym(:,1)
      pcesym => spcsym(:,2)
      fixsym => spcsym(:,3)
      grpsym => spcsym(:,4)
      colsym => spcsym(:,5)
      new_solsym => spcsym(:,6)
      grp_mem_sym => spcsym(:,7)

      slvdsym => spcsym(:,8)

      end subroutine VAR_INI

      end module VAR_MOD

      module RXT_MOD
!-----------------------------------------------------------------------
!        ... Mozart reaction variables
!-----------------------------------------------------------------------

      implicit none

      integer :: rxt_lim
      integer :: rxtnt_lim
      integer :: prd_lim, prd_limp1
      integer ::   phtcnt = 0, &                    ! count of photolysis reactions
                   hetcnt = 0, &                    ! count of heterogeneous processes
                   usrcnt = 0, &                    ! count of "extraneous" forcing processes
                   rxntot = 0, &                    ! count of photo and gas phase reactions
                   gascnt = 0                       ! count of gas phase reactions
      integer, allocatable :: &
                   fixmap(:,:,:), &
                   prdmap(:,:)
      integer, dimension(2) :: &
                   fixcnt = 0, &                    ! count of reactions with fixed rxtnts
                   rxmcnt = 0, &                    ! count of reactions with sol rxtnts
                   ipcel = 0                        ! not used
      integer, dimension(3) :: &
                   ipcep  = 0                       ! not used
      integer  ::  prdcnt = 0, &                    ! entries in prdmap matrix
                   rxpcnt = 0, &                    ! entries in rxparm matrix
                   troecnt = 0                      ! count of troe rates

      integer, allocatable :: &
                   rxmap(:,:,:), &
                   pcel(:,:,:), &                   ! not used
                   pcep(:,:,:)                      ! not used

      integer, allocatable :: &
                   rxptab(:), &
                   troetab(:), &
                   hetmap(:,:), &
                   usrmap(:)

      integer, dimension(2) :: &
                   grp_rat_cnt = 0, &
                   rel_rxt_cnt = 0
      integer, allocatable :: &
                   grp_rat_map(:,:,:), &
                   rxt_to_grp_map(:,:), &
                   rel_rxt_pntr(:,:), &
                   rel_rxt_map(:,:,:)

      integer  ::  pcoeff_cnt = 0
      integer, allocatable :: &
                   pcoeff_ind(:)
      real, allocatable :: &
                   pcoeff(:,:), &
                   rxparm(:,:), &
                   troe_rxparm(:,:)
      character(len=16), allocatable :: &
                   sym_rates(:,:), &
                   troe_sym_rates(:,:), &
                   pht_alias(:,:), &
                   pht_alias_mult(:,:), &
                   rxt_tag(:)
      character(len=16), allocatable :: &
                   phtsym(:)
      logical, allocatable :: &
                   rxt_has_tag(:), &
                   rxt_has_alias(:) 
      logical, allocatable :: &
                   cph_flg(:) 
      logical, allocatable :: &
                   frc_from_dataset(:)

      integer :: cls_rxt_cnt(4,5) = 0
      integer, allocatable :: &
                   cls_rxt_map(:,:,:)

      contains

      subroutine RXT_INI()
!-----------------------------------------------------------------------
!        ... Allocate and initialize reaction variables
!-----------------------------------------------------------------------

      use VAR_MOD, only : var_lim

      implicit none

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer :: astat

!-----------------------------------------------------------------------
!        ... Set reaction limits
!-----------------------------------------------------------------------
      rxt_lim   = 900
      rxtnt_lim = 3
!     prd_lim   = 16
      prd_lim   = 24
      prd_limp1 = prd_lim + 1

      allocate( fixmap(var_lim,3,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate fixmap'
	 stop
      end if
      fixmap(:,:,:) = 0
      allocate( prdmap(var_lim,prd_limp1),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate prdmap'
	 stop
      end if
      prdmap(:,:) = 0
      allocate( rxmap(rxt_lim,prd_lim+3,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate rxmap'
	 stop
      end if
      rxmap(:,:,:) = 0
      allocate( pcel(var_lim,6,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate pcel'
	 stop
      end if
      pcel(:,:,:) = 0
      allocate( pcep(var_lim,5,3),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate pcep'
	 stop
      end if
      pcep(:,:,:) = 0
      allocate( rxptab(rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate rxptab'
	 stop
      end if
      rxptab(:) = 0
      allocate( troetab(rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate troetab'
	 stop
      end if
      troetab(:) = 0
      allocate( hetmap(rxt_lim,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate hetmap'
	 stop
      end if
      hetmap(:,:) = 0
      allocate( usrmap(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate usrmap'
	 stop
      end if
      usrmap(:) = 0
      allocate( grp_rat_map(rxt_lim,3,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate grp_rat_map'
	 stop
      end if
      grp_rat_map(:,:,:) = 0
      allocate( rxt_to_grp_map(rxt_lim,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate rxt_to_grp_map'
	 stop
      end if
      rxt_to_grp_map(:,:) = 0
      allocate( rel_rxt_pntr(rxt_lim,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate rel_rxt_pntr'
	 stop
      end if
      rel_rxt_pntr(:,:) = 0
      allocate( rel_rxt_map(rxt_lim,3,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate rel_rxt_map'
	 stop
      end if
      rel_rxt_map(:,:,:) = 0
      allocate( pcoeff_ind(rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate pcoeff_ind'
	 stop
      end if
      pcoeff_ind(:) = 0
      allocate( pcoeff(prd_lim,rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate pcoeff'
	 stop
      end if
      pcoeff(:,:) = 0.
      allocate( rxparm(2,rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate rxparm'
	 stop
      end if
      rxparm(:,:) = 0.
      allocate( troe_rxparm(5,rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate troe_rxparm'
	 stop
      end if
      troe_rxparm(:,:) = 0.
      allocate( sym_rates(2,rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate sym_rates'
	 stop
      end if
      sym_rates(:,:) = ' '
      allocate( troe_sym_rates(5,rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate troe_sym_rates'
	 stop
      end if
      troe_sym_rates(:,:) = ' '
      allocate( phtsym(rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate phtsym'
	 stop
      end if
      phtsym(:) = ' '
      allocate( rxt_tag(rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate rxt_tag'
	 stop
      end if
      rxt_tag(:) = ' '
      allocate( rxt_has_tag(rxt_lim),rxt_has_alias(rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate rxt_has_tag,rxt_has_alias'
	 stop
      end if
      rxt_has_tag(:)   = .false.
      rxt_has_alias(:) = .false.
      allocate( cph_flg(rxt_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate cph_flg'
	 stop
      end if
      cph_flg(:) = .false.
      allocate( pht_alias(rxt_lim,2),pht_alias_mult(rxt_lim,2),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate pht_alias,pht_alias_mult'
	 stop
      end if
      pht_alias(:,:) = ' '
      pht_alias_mult(:,:) = '1.'
      allocate( frc_from_dataset(var_lim),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate frc_from_dataset'
	 stop
      end if
      frc_from_dataset(:) = .false.
      allocate( cls_rxt_map(rxt_lim,prd_lim+3,5),stat=astat )
      if( astat /= 0 ) then
	 write(*,*) 'RXT_INI: Failed to allocate cls_rxt_map'
	 stop
      end if
      cls_rxt_map(:,:,:) = 0

      end subroutine RXT_INI

      end module RXT_MOD
