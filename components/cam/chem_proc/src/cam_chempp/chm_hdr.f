
      subroutine chm_hdr( rxt_tag_cnt, hetcnt, usrcnt, cls_rxt_cnt, radj_flag, phtcnt, &
                          rxpcnt, rxparm, rxntot, ncol, nfs, nslvd, &
                          indexm, indexh2o, spcno, relcnt, grpcnt, &
                          clscnt, iter_counts, nzcnt, vec_ftns, machine, chemistry )
!-----------------------------------------------------------------------
!        ... Write the chemistry "header" file
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::    rxt_tag_cnt
      integer, intent(in) ::    hetcnt                ! count of washout processes
      integer, intent(in) ::    usrcnt                ! count of extraneous forcing
      integer, intent(in) ::    phtcnt                ! count of photorates
      integer, intent(in) ::    rxpcnt                ! count of specified rates
      integer, intent(in) ::    rxntot                ! count of totol reactions
      integer, intent(in) ::    ncol                  ! number of column integrals
      integer, intent(in) ::    nfs                   ! number of "fixed" species
      integer, intent(in) ::    nslvd                 ! number of "short lived" species
      integer, intent(in) ::    indexm                ! index for "m"
      integer, intent(in) ::    indexh2o              ! index for h2o
      integer, intent(in) ::    spcno                 ! total number of xported species
      integer, intent(in) ::    relcnt                ! number of "relative" species
      integer, intent(in) ::    grpcnt                ! number of group species
      integer, intent(in) ::    nzcnt(2)              ! number of non-zero entries in lu
      integer, intent(in) ::    clscnt(5)             ! solution class count
      integer, intent(in) ::    iter_counts(4)        ! iteration counts
      integer, intent(in) ::    cls_rxt_cnt(4,5)      ! class reaction count

      real, intent(in)    ::    rxparm(2,*)           ! rxtn rate parms
      logical, intent(in) ::    radj_flag             ! rxt adjust flag
      logical, intent(in) ::    vec_ftns              ! vector function flag
      logical, intent(in) ::    chemistry             ! chemistry flag

      character(len=16), intent(in)  :: machine        ! target machine

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer  ::    gascnt              ! number of gas phase rxtns
      logical  ::  lexist

      inquire( file = 'chem.h', exist = lexist )
      if( lexist ) then
	 call system( 'rm chem.h' )
      end if
      open( unit = 30, file = 'chem.h' )

      write(30,'(''# define RXTTAGCNT '',i5)') rxt_tag_cnt
      write(30,'(''# define HETCNT '',i5)') hetcnt
      write(30,'(''# define EXTCNT '',i5)') usrcnt
      gascnt = sum( cls_rxt_cnt(1,1:5) )
      write(30,'(''# define CLSINDPRD '',i5)') gascnt
      write(30,'(''# define CLSINDPRD1 '',i5)') cls_rxt_cnt(1,1)
      write(30,'(''# define CLSINDPRD2 '',i5)') cls_rxt_cnt(1,2)
      write(30,'(''# define CLSINDPRD3 '',i5)') cls_rxt_cnt(1,3)
      write(30,'(''# define CLSINDPRD4 '',i5)') cls_rxt_cnt(1,4)
      write(30,'(''# define CLSINDPRD5 '',i5)') cls_rxt_cnt(1,5)
      write(30,'(''# define IMP_NZCNT '',i5)') nzcnt(1)
      write(30,'(''# define ROD_NZCNT '',i5)') nzcnt(2)
      gascnt = cls_rxt_cnt(2,4) + cls_rxt_cnt(4,4)
      write(30,'(''# define IMP_LINCNT '',i5)') gascnt
      gascnt = cls_rxt_cnt(2,5) + cls_rxt_cnt(4,5)
      write(30,'(''# define ROD_LINCNT '',i5)') gascnt
      write(30,'(''# define IMP_NLNCNT '',i5)') cls_rxt_cnt(3,4)
      write(30,'(''# define ROD_NLNCNT '',i5)') cls_rxt_cnt(3,5)
      if( radj_flag ) then
         write(30,'(''# define RADJFLAG'')')
      end if
      write(30,'(''# define PHTCNT '',i5)') phtcnt
      write(30,'(''# define PHTCNTP1 '',i5)') phtcnt+1
      write(30,'(''# define RXNCNT '',i5)') rxntot
      gascnt = rxntot - phtcnt
      write(30,'(''# define GASCNT '',i5)') gascnt
      write(30,'(''# define SETRXNCNT '',i5)') rxpcnt
      write(30,'(''# define USRRXNCNT '',i5)') gascnt - rxpcnt
      gascnt = count( rxparm(2,1:rxpcnt) /= 0. )
      write(30,'(''# define TDEPCNT '',i5)') gascnt
      write(30,'(''# define NCOL '',i5)') ncol
      write(30,'(''# define NFS '',i5)') nfs
      write(30,'(''# define NSLVD '',i5)') nslvd
      write(30,'(''# define INDEXM '',i5)') indexm
      write(30,'(''# define INDEXH2O '',i5)') indexh2o
      write(30,'(''# define PCNST '',i5)') spcno
      write(30,'(''# define PCNSTP2 '',i5)') spcno+2
      write(30,'(''# define RELCNT '',i5)') relcnt
      write(30,'(''# define GRPCNT '',i5)') grpcnt
      write(30,'(''# define CLSCNT1 '',i5)') clscnt(1)
      write(30,'(''# define CLSCNT2 '',i5)') clscnt(2)
      write(30,'(''# define CLSCNT3 '',i5)') clscnt(3)
      write(30,'(''# define CLSCNT4 '',i5)') clscnt(4)
      write(30,'(''# define CLSCNT5 '',i5)') clscnt(5)
      write(30,'(''# define EBIITERMAX '',i5)') MAX( 1,iter_counts(4) )
      write(30,'(''# define HOVITERMAX '',i5)') MAX( 1,iter_counts(1) )
      write(30,'(''# define IMPITERMAX '',i5)') MAX( 1,iter_counts(2) )
      write(30,'(''# define IMPJACITER '',i5)') MAX( 1,iter_counts(3) )
      if( chemistry ) then
         write(30,'(''# define TROP_CHEM'')')
      end if
      close( 30 )
      
      end subroutine chm_hdr
