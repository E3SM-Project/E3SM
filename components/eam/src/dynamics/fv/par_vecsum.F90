!-----------------------------------------------------------------------
!BOP
! !ROUTINE: par_vecsum --- Calculate vector sum bit-wise consistently
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine par_vecsum(jm, jfirst, jlast, InVector, te0,    &
                            incomm, npryuse)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
#if defined ( SPMD )
      use parutilitiesmodule, only : parexchangevector
#endif
      use shr_kind_mod, only: r8 => shr_kind_r8

      implicit none

! !INPUT PARAMETERS:
      integer jm                    ! global latitudes
      integer jfirst                ! first latitude on this PE
      integer jlast                 ! last latitude on this PE
      real (r8) InVector(jm)        ! input vector to be summed
      integer incomm                ! communicator
      integer npryuse               ! number of subdomains

! !OUTPUT PARAMETERS:
      real (r8) te0                 ! sum of all vector entries

! !DESCRIPTION:
!     This subroutine calculates the sum of InVector in a reproducible
!     (sequentialized) fashion which should give bit-wise identical
!     results irrespective of the number of MPI processes.
!
! !CALLED FROM:
!     te_map and benergy
!
! !REVISION HISTORY:
!
!     BWS 00.01.15 : Created
!     WS  00.06.02 : Replaced MPI calls with ParExchangeVector; docu.
!     WS  00.08.29 : SPMD instead of MPI_ON
!     AM  01.06.15 : general communicator
!
!EOP
!---------------------------------------------------------------------
!BOC
 
! !Local
      real(r8), parameter ::  D0_0                    =  0.0_r8
      real (r8) tte_all(jm)
      integer j

#if defined ( SPMD )
      real (r8)     tte_send(npryuse*(jlast-jfirst+1))
      integer  sendcount(npryuse)
      integer  recvcount(npryuse)
      integer  ipe, icount
#endif

      te0 = D0_0
#if defined ( SPMD ) 
      icount=0
      do ipe=1,npryuse
        sendcount(ipe) = jlast-jfirst+1
        do j=jfirst, jlast
          icount=icount+1    
          tte_send(icount)=InVector(j)
        enddo
      enddo
      call parexchangevector( incomm, sendcount, tte_send,          &
                              recvcount, tte_all )
#else
      do j=1, jm
        tte_all(j)=InVector(j)
      enddo
#endif

      te0 = D0_0
      te0 = te0 + tte_all(1)     !in oder to compare to SMP-only
      te0 = te0 + tte_all(jm)    !in oder to compare to SMP-only

      do j=2,jm-1
        te0 = te0 + tte_all(j)
      enddo

      return
!EOC
      end
!-----------------------------------------------------------------------

