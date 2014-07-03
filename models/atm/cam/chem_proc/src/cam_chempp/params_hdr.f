
      subroutine PARAMS_HDR( plon, plonl, plat, plev, &
			     phtcnt, rxntot, &
			     adv_cnt, nadv_cnt, histout_cnt, &
		             chemistry, diffusion, convection, arch_type, &
                             filespec )
!-----------------------------------------------------------------------
!        ... Make the params.h file
!-----------------------------------------------------------------------

      use IO, only : lout

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::     plon, plonl, plat, plev
      integer, intent(in) ::     phtcnt, rxntot
      integer, intent(in) ::     adv_cnt, nadv_cnt
      integer, intent(in) ::     histout_cnt(20,2)
      logical, intent(in) ::     chemistry, diffusion, convection
      character(len=16), intent(in) :: arch_type
      character(len=*), intent(in) :: filespec

!-----------------------------------------------------------------------
!        ... The local variables
!-----------------------------------------------------------------------
      integer           :: ios
      logical           :: lexist
      
      INQUIRE( file = TRIM(filespec), exist = lexist )
      if( lexist ) then
	 call SYSTEM( 'rm '//TRIM(filespec) )
      end if
      OPEN( unit = 30, file = TRIM(filespec), iostat=ios )
      if( ios /= 0 ) then
	 write(lout,*) 'PARAMS_HDR: Failed to open simulation datafile ',TRIM(filespec),' ;error = ',ios
	 stop
      end if
      write(30,'(a)') '# ifndef PARAMS_H'
      write(30,'(a)') '# define PARAMS_H'
      write(30,'('' '')')
      write(30,'(a)') '# define CALC_ETADOT'
      if( diffusion ) then
         write(30,'(a)') '# define DI_VDIFF'
      else
         write(30,'(a)') '# define AR_VDIFF'
      end if
      if( convection ) then
         write(30,'(a)') '# define DI_CONV_CCM'
      else
         write(30,'(a)') '# define AR_CONV_CCM'
      end if
      write(30,'(a)') '# define DI_CLOUD_PHYS'
      if( chemistry ) then
         write(30,'(a)') '# define TROP_CHEM'
      end if
      write(30,'('' '')')
      write(30,'(a)') '# endif'
      
      end subroutine PARAMS_HDR
