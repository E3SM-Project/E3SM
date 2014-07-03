
      subroutine RES_HDR( plon, &
                          plonl, &
                          plat, &
                          plev, &
                          jintmx, &
                          nxpt, &
                          arch_type, &
                          cpucnt )

      implicit none

!-----------------------------------------------------------------------
!        ... The arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::    plon
      integer, intent(in) ::    plonl
      integer, intent(in) ::    plat
      integer, intent(in) ::    plev
      integer, intent(in) ::    jintmx
      integer, intent(in) ::    nxpt
      integer, intent(in) ::    cpucnt
      character(len=16), intent(in) :: arch_type

!-----------------------------------------------------------------------
!        ... The local variables
!-----------------------------------------------------------------------
      character(len=72) :: comment
      logical      :: lexist
      
      INQUIRE( file = 'res.h', exist = lexist )
      if( lexist ) then
	 call SYSTEM( 'rm res.h' )
      end if
      OPEN( unit = 30, file = 'res.h' )

      write(30,'(''# define PLON '',i3)') plon
      write(30,'(''# define PLONP2 '',i3)') plon + 2
      if( arch_type == 'HYBRID' ) then
         write(30,'(''# define PLONL '',i3)') plonl
      end if
      write(30,'(''# define NXPT '',i3)') nxpt
      write(30,'(''# define JINTMX '',i3)') jintmx
      write(30,'(''# define NXPTJ '',i3)') nxpt + jintmx
      write(30,'(''# define PLOND '',i3)') plon + 1 + 2*nxpt
      write(30,'(''# define PLAT '',i3)') plat
      write(30,'(''# define PLEV '',i3)') plev
      write(30,'(''# define PLEVP '',i3)') plev + 1
      write(30,'(''# define PLEVM '',i3)') plev - 1
      write(30,'(''# define PLNPLV '',i6)') plon*plev
      write(30,'(''# define I1 '',i3)') 1 + nxpt
      write(30,'(''# define I1M '',i3)') nxpt
      write(30,'(''# define J1 '',i3)') 1 + nxpt + jintmx
      write(30,'(''# define J1M '',i3)') nxpt + jintmx
      write(30,'(''# define PTIML 2 '')')
      
      CLOSE(30)
      
      end subroutine RES_HDR
