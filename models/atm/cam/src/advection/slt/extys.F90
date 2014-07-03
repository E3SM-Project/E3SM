
subroutine extys(pkcnst  ,pkdim   ,fb      ,kloop)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Fill latitude extensions of a scalar extended array and
! Copy data to the longitude extensions of the extended array
! 
! Method: 
! This is done in 2 steps:
!   1) interpolate to the pole points; use the mean field value on the
!      Gaussian latitude closest to the pole.
!   2) add latitude lines beyond the poles.
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plat
  use scanslt,      only: nxpt, plond, beglatex, endlatex, platd, nlonex, &
                          jintmx
  implicit none

!------------------------------Parameters-------------------------------
  integer, parameter :: istart = nxpt+1           ! index to start computation
  integer, parameter :: js = 1    + nxpt + jintmx ! index of southernmost model lat
  integer, parameter :: jn = plat + nxpt + jintmx ! index of northernmost model lat
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
  integer , intent(in) :: pkcnst   ! dimensioning construct for 3-D arrays
  integer , intent(in) :: pkdim    ! vertical dimension
  real(r8), intent(inout) :: fb(plond,pkdim*pkcnst,beglatex:endlatex) ! Output is same as on entry 
                          !except with the pole latitude and extensions beyond it filled.
  integer, intent(in) :: kloop ! If you want to limit the extent of looping over pcnst
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,j,k             ! indices
  integer istop             ! index to stop  computation
  integer nlon2             ! half the number of real longitudes
  real(r8) zave             ! accumulator for zonal averaging
  integer pk                ! dimension to loop over
!-----------------------------------------------------------------------
!
! Fill north pole line.
!
  pk = pkdim*kloop
#if ( defined SPMD )
  if (jn+1<=endlatex) then
#endif
     do k=1,pkdim*pkcnst
        zave = 0.0_r8
        istop  = nxpt + nlonex(jn)
        do i=istart,istop
           zave = zave + fb(i,k,jn  )
        end do
        zave = zave/nlonex(jn)
        istop  = nxpt + nlonex(jn+1)
        do i=istart,istop
           fb(i,k,jn+1) = zave
        end do
     end do
#if ( defined SPMD )
  end if
#endif
!
! Fill northern lines beyond pole line.
!
  if( jn+2 <= platd )then
     do j=jn+2,platd
#if ( defined SPMD )
        if (j<=endlatex) then
#endif
           nlon2 = nlonex(j)/2
           do k=1,pk
              do i=istart,istart+nlon2-1
                 fb(      i,k,j) = fb(nlon2+i,k,2*jn+2-j)
                 fb(nlon2+i,k,j) = fb(      i,k,2*jn+2-j)
              end do
           end do
#if ( defined SPMD )
        end if
#endif
     end do
  end if
!
! Fill south pole line.
!
#if ( defined SPMD )
  if (js-1>=beglatex) then
#endif
     do k=1,pk
        zave = 0.0_r8
        istop  = nxpt + nlonex(js)
        do i = istart,istop
           zave = zave + fb(i,k,js  )
        end do
        zave = zave/nlonex(js)
        istop  = nxpt + nlonex(js-1)
        do i=istart,istop
           fb(i,k,js-1) = zave
        end do
     end do
#if ( defined SPMD )
  end if
#endif
!
! Fill southern lines beyond pole line.
!
  if( js-2 >= 1 )then
     do j=1,js-2
#if ( defined SPMD )
        if (j>=beglatex) then
#endif
           nlon2 = nlonex(j)/2
           do k=1,pk
              do i=istart,istart+nlon2-1
                 fb(      i,k,j) = fb(nlon2+i,k,2*js-2-j)
                 fb(nlon2+i,k,j) = fb(      i,k,2*js-2-j)
              end do
           end do
#if ( defined SPMD )
        end if
#endif
     end do
  end if

  return
end subroutine extys
