
      module mo_ver_hdr

      contains

      subroutine ver_hdr( options, &
                          plon, plonl, plev, &
                          machine, &
                          model, &
                          arch_type, &
                          ohstflag, &
                          diagprnt, &
                          tavgprnt, &
                          srf_flx_cnt, &
		          hetcnt, rxntot, clscnt, nzcnt, spcno, &
                          dvel_cnt )

      implicit none

!-----------------------------------------------------------------------
!        ... The arguments
!-----------------------------------------------------------------------
      integer, intent(in)  ::      srf_flx_cnt      ! species with srf flux
      integer, intent(in)  ::      dvel_cnt         ! species with dep vel
      integer, intent(in)  ::      plon, plonl, plev
      integer, intent(in)  ::      hetcnt, rxntot, spcno
      integer, intent(in)  ::      nzcnt(2)
      integer, intent(in)  ::      clscnt(5)
      character(len=16), intent(in) :: machine       ! target machine
      character(len=16), intent(in) :: model         ! target model
      character(len=16), intent(in) :: arch_type     ! architecture
      logical, intent(in)  ::      options(*)       ! options array
      logical, intent(in)  ::      ohstflag         ! hist tape write flag
      logical, intent(in)  ::      diagprnt         ! chktrc, negtrc diag printout flag
      logical, intent(in)  ::      tavgprnt         ! time averaged printout flag

!-----------------------------------------------------------------------
!        ... The local variables
!-----------------------------------------------------------------------
      integer           :: i, cache_factor
      integer           :: up_bound(2)
      logical           :: lexist
      
      inquire( file = 'version.h', exist = lexist )
      if( lexist ) then
	 call system( 'rm version.h' )
      end if
      open( unit = 30, file = 'version.h' )

      if( options(1) ) then
         write(30,'(''# define CHEM'')')
      end if

      if( options(2) ) then
         write(30,'(''# define CRAY'')')
      end if

      if( ohstflag ) then
         write(30,'(''# define HISTTAPE'')')
      end if

      if( diagprnt ) then
         write(30,'(''# define DIAGPRNT'')')
      end if

      if( tavgprnt ) then
         write(30,'(''# define TAVGPRNT'')')
      end if

      if( options(12) ) then
         write(30,'(''# define RXTNLOOKUP'')')
      end if

      if( options(14) ) then
         write(30,'(''# define F90'')')
      end if

      if( options(16) ) then
         write(30,'(''# define USRHOOK'')')
      end if

      if( options(17) ) then
         write(30,'(''# define MODULES'')')
      end if

      if( srf_flx_cnt /= 0 ) then
         write(30,'(''# define SFLUX'')')
      end if

      if( dvel_cnt /= 0 ) then
         write(30,'(''# define DVEL'')')
      end if

      select case( machine )
         case( 'INTEL' )
            write(30,'(''# define CLSZE 1'')')
            write(30,'(''# define MACHINE_INTEL'')')
         case( 'ALPHA', 'IBM' )
	    do i = 1,2
	       up_bound(i) = 2*nzcnt(i) + rxntot + hetcnt + 6*clscnt(i+3) + spcno
	       up_bound(i) = CEILING( 8.*1024./REAL(up_bound(i)) )
	    end do
	    i = MINVAL( up_bound(:) )
	    up_bound(1) = i
	    do i = 2,up_bound(1)
	       if( MOD( plonl,i) == 0 ) then
		  cache_factor = i
	       end if
	    end do
            write(30,'(''# define CLSZE '',i3)')  MAX( 1,cache_factor )
	    if( machine == 'IBM' ) then
               write(30,'(''# define MACHINE_IBM'')')
	    else
               write(30,'(''# define MACHINE_ALPHA'')')
	    end if
         case( 'CRAY': 'CRAYYMP', 'J90', 'C90' )
            write(30,'(''# define CLSZE '',i3)') plon
            write(30,'(''# define MACHINE CRAY'')')
         case( 'NEC', 'FUJITSU' )
            write(30,'(''# define CLSZE '',i5)') plon*plev
	    if( machine == 'NEC' ) then
               write(30,'(''# define MACHINE_NEC'')')
	    else
               write(30,'(''# define MACHINE_FUJITSU'')')
	    end if
         case default
	    if( arch_type == 'HYBRID' ) then
               write(30,'(''# define CLSZE  4'')')
	    else
               write(30,'(''# define CLSZE '',i3)') plon
            end if
      end select
      if( model == 'MOZART' ) then
         write(30,'(''# define MOZART '')')
      else if( model == 'MOZART' ) then
         write(30,'(''# define CAM '')')
      end if

      close(30)
      
      end subroutine ver_hdr

      end module mo_ver_hdr
