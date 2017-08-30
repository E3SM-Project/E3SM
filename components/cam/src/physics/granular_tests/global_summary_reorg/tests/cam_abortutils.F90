! cam_abortutils.F90
!    This F90 module file is a special version of the equivalent ACME (and CAM5) module.
!    It provides the functionality needed by the cambox offline code
!    that is used for development and testing of the modal aerosol module (MAM),
!    but (in most cases) not all the functionality of the equivalent ACME module.
!    Also, it may have been taken from a version of CAM5 that was older
!    than ACME-V0 (i.e., pre 2014).

      module cam_abortutils

      use cam_logfile, only:  iulog

      implicit none

      public

      contains

      subroutine endrun( msg )
      character(len=*), optional, intent(in) :: msg
      integer :: lunout
      lunout = 6
      if ( present( msg ) ) then
         if (msg /= ' ') then
            write(lunout,'(/a/a/)') &
               '*** stopping in ENDRUN with msg =', trim(msg)
         else
            write(lunout,'(/a/)') &
               '*** stopping in ENDRUN with msg = blank'
         end if
      else
         write(lunout,'(/a/)') &
            '*** stopping in ENDRUN with no msg'
      end if
      stop
      end subroutine endrun

      end module cam_abortutils
