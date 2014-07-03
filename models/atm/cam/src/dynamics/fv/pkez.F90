!-----------------------------------------------------------------------
!BOP
! !ROUTINE: pkez --- Calculate solution to hydrostatic equation
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine pkez(nx, im, km, jfirst, jlast, kfirst, klast,    &
                      ifirst, ilast, pe, pk, akap, ks, peln, pkz, eta)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8

      implicit none

!
! This routine may be called assuming either yz or xy decompositions.
! For xy decomposition, the effective "nx" is 1.
!

! !INPUT PARAMETERS:
      integer nx                          ! SMP decomposition in x
      integer im, km                      ! Dimensions
      integer jfirst, jlast               ! Latitude strip
      integer kfirst, klast               ! Vertical strip
      integer ifirst, ilast               ! Longitude strip
      real (r8)  pe(ifirst:ilast, kfirst:klast+1, jfirst:jlast)    ! Edge pressure
      integer ks
      logical eta     ! Is on ETA coordinate?
                      ! True:  input pe    ; output pk, pkz, peln
                      ! False: input pe, pk; output     pkz, peln
      real (r8) akap

! !INPUT/OUTPUT PARAMETERS:
      real (r8)  pk(ifirst:ilast,jfirst:jlast,kfirst:klast+1)

! !OUTPUT PARAMETERS
      real (r8) pkz(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real (r8) peln(ifirst:ilast, kfirst:klast+1, jfirst:jlast)   ! log pressure (pe) at layer edges

! !DESCRIPTION:
!
!
! !CALLED FROM:
!     te_map and fvccm3
!
! !REVISION HISTORY:
!
!     WS  99.05.19 : Removed fvcore.h
!     WS  99.07.27 : Limited region to jfirst:jlast
!     WS  99.10.22 : Deleted cp as argument (was not used)
!     WS  99.11.05 : Documentation; pruning of arguments
!     SJL 00.01.02 : SMP decomposition in i
!     AAM 00.08.10 : Add kfirst:klast
!     AAM 01.06.27 : Add ifirst:ilast
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local
      real (r8) pk2(ifirst:ilast, kfirst:klast+1)
      real (r8) pek
      real (r8) lnp
      integer i, j, k, itot, nxu
      integer ixj, jp, it, i1, i2

      itot = ilast - ifirst + 1
! Use smaller block sizes only if operating on full i domain
      nxu = 1
      if (itot .eq. im) nxu = nx

      it = itot / nxu
      jp = nxu * ( jlast - jfirst + 1 )

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(ixj, i1, i2, i, j, k, pek, lnp, pk2)

! WS 99.07.27 : Limited region to jfirst:jlast

      do 1000 ixj=1,jp

         j  = jfirst + (ixj-1) / nxu
         i1 = ifirst + it * mod(ixj-1, nxu)
         i2 = i1 + it - 1

        if ( eta ) then

! <<<<<<<<<<< Eta cordinate Coordinate  >>>>>>>>>>>>>>>>>>>
          if (kfirst .eq. 1) then
            pek =     pe(i1,1,j)**akap
            lnp = log(pe(i1,1,j))

            do i=i1,i2
               pk2(i,1)   = pek
              peln(i,1,j) = lnp
            enddo
          endif

          if(ks .ne. 0) then
            do k=max(2,kfirst), min(ks+1,klast+1)
              pek = pe(i1,k,j)**akap
              lnp = log(pe(i1,k,j))
              do i=i1,i2
                 pk2(i,k)   = pek
                peln(i,k,j) =  lnp
              enddo
            enddo

            do k=kfirst, min(ks,klast)
              pek = (       pk2(i1,k+1)   - pk2(i1,k))   /     &
                    (akap*(peln(i1,k+1,j) - peln(i1,k,j)) )
              do i=i1,i2
                 pkz(i,j,k) = pek
              enddo
            enddo
          endif

          do k=max(ks+2,kfirst), klast+1
#if !defined( VECTOR_MATH )
            do i=i1,i2
               pk2(i,k) = pe(i,k,j)**akap
            enddo
#else
            call vlog(pk2(i1,k), pe(i1,k,j), it)
            do i=i1,i2
               pk2(i,k) = akap * pk2(i,k)
            enddo
            call vexp(pk2(i1,k), pk2(i1,k), it)
#endif
          enddo

          do k=max(ks+2,kfirst), klast+1
            do i=i1,i2
               peln(i,k,j) =  log(pe(i,k,j))
            enddo
          enddo

          do k=max(ks+1,kfirst), klast
            do i=i1,i2
               pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k)) /         &
                            (akap*(peln(i,k+1,j) - peln(i,k,j)) )
            enddo
          enddo

          do k=kfirst, klast+1
            do i=i1,i2
               pk(i,j,k) = pk2(i,k)
            enddo
          enddo

        else

! <<<<<<<<<<< General Coordinate  >>>>>>>>>>>>>>>>>>>

          if (kfirst .eq. 1) then
            pek =     pk(i1,j,1)
            lnp = log(pe(i1,1,j))

            do i=i1,i2
               pk2(i,1) = pek
               peln(i,1,j) = lnp
            enddo
          endif

          do k=max(2,kfirst), klast+1
             do i=i1,i2
                peln(i,k,j) =  log(pe(i,k,j))
                pk2(i,k) =  pk(i,j,k)
             enddo
          enddo

          do k=kfirst, klast
             do i=i1,i2
                pkz(i,j,k) = (       pk2(i,k+1) - pk2(i,k) )  /    &
                             (akap*(peln(i,k+1,j) - peln(i,k,j)) )
             enddo
          enddo

        endif
1000  continue

      return
!EOC
      end
!-----------------------------------------------------------------------
