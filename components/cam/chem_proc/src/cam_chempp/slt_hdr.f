      subroutine SLT_HDR( cray, &
                          multitask, &
                          cpucnt, &
                          machine )

      implicit none

!-----------------------------------------------------------------------
!        ... The arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::     cpucnt
      logical, intent(in) ::     cray
      logical, intent(in) ::     multitask
      character(len=16), intent(in) :: machine

!-----------------------------------------------------------------------
!        ... The local variables
!-----------------------------------------------------------------------
      logical           :: lexist
      
      INQUIRE( file = 'slt.h', exist = lexist )
      if( lexist ) then
	 call SYSTEM( 'rm slt.h' )
      end if
      OPEN( unit = 30, file = 'slt.h' )

      if( cray .and. multitask ) then
         write(30,'(''# define MT'')')
      else if( .not. cray ) then
         write(30,'(''# define NOCRAY'')')
         write(30,'(''# define PORT'')')
	 if( multitask ) then
            write(30,'(''# define MPP'')')
	    write(30,'(''# define NCPUS '',i3)') cpucnt
	 end if
      end if
      write(30,'(''# define '',a8)') machine

      CLOSE(30)
      
      end subroutine SLT_HDR
