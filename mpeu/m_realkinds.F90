!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_realkinds - real KIND definitions
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_realkinds
      implicit none
      private	! except

      public :: kind_r4		! real*4
      public :: kind_r8		! real*8
      public :: kind_r		! default real
      public :: SP		! default REAL
      public :: DP		! default DOUBLE_PRECISION

      real*4,parameter :: R4=1.
      real*8,parameter :: R8=1.
      real,  parameter :: R =1.

      integer,parameter :: SP = kind(1.  )
      integer,parameter :: DP = kind(1.D0)

      integer,parameter :: kind_r4=kind(R4)
      integer,parameter :: kind_r8=kind(R8)
      integer,parameter :: kind_r =kind(R )

! !REVISION HISTORY:
! 	19Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP
!_______________________________________________________________________
  character(len=*),parameter :: myname='m_realkinds'

end module m_realkinds
