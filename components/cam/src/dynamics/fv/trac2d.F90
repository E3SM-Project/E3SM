!-----------------------------------------------------------------------
!BOP
! !ROUTINE: trac2d --- Remap Lagrangian to fixed coordinates
!
! !INTERFACE:
 subroutine trac2d( grid,    dp1,  tracer, cx,     cy,      &
                    mfx,     mfy,   iord,  jord,   fill,    &
                    nlo,     nhi,   va,    flx )
 
! !USES:

   use shr_kind_mod, only  : r8 => shr_kind_r8, r4 => shr_kind_r4
   use tp_core, only       : tp2c
   use fill_module, only   : fillxy
   use dynamics_vars, only : T_FVDYCORE_GRID
   use FVperf_module, only : FVstartclock, FVstopclock, FVbarrierclock

#if defined( SPMD )
   use parutilitiesmodule, only: maxop, parcollective
   use mod_comm, only : mp_send4d_ns, mp_recv4d_ns,  &
                        mp_send4d_ns_r4, mp_recv4d_ns_r4,        &
                        mp_send3d_2, mp_recv3d_2
#endif

   implicit none

! !INPUT PARAMETERS:

   type (T_FVDYCORE_GRID), intent(inout) :: grid
   integer, intent(in):: iord,  jord

   logical, intent(in):: fill
   integer, intent(in):: nlo,  nhi   ! Tracer index range

! !INPUT/OUTPUT PARAMETERS:
   real(r8), intent(inout):: dp1(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
   real(r8), intent(inout):: cx(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d,grid%kfirst:grid%klast)
   real(r8), intent(inout):: cy(grid%im,grid%jfirst:grid%jlast+1,grid%kfirst:grid%klast)
   real(r8), intent(inout):: mfx(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
   real(r8), intent(inout):: mfy(grid%im,grid%jfirst:grid%jlast+1,grid%kfirst:grid%klast)
   real(r8), intent(inout):: tracer(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast,grid%ntotq)

! !OUTPUT PARAMETERS:
   real(r8), intent(out):: va(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
   real(r8), intent(out):: flx(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)

! !DESCRIPTION:
!
!  Perform large-time-step tracer transport using accumulated Courant
!  numbers (cx, cy) and the mass fluxes (mfx, mfy) within the Lagrangian
!  layers.  This routine is 100\% parallel in the vertical direction
!  (with SMP).  Merdional Courant number will be further split, if
!  necessary, to ensure stability.  Cy <= 1 away from poles; Cy $\le$
!  1/2 at the latitudes closest to the poles.
!
! !REVISION HISTORY:
!
!   SJL 99.04.13:  Delivery
!   WS  99.05.26:  Added jfirst:jlast concept; im, jm, km as parameters
!                  replaced IMR, JMR, JNP, NL with IM, JM-1, JM and KM
!   WS  99.09.27:  Documentation; indentation; jfirst:jlast 
!   WS  99.09.30:  Ghosting; loop limits; full parallelization; tested
!   SJL 99.10.15:  nsplt migrated to outermost loop to remove bug
!   SJL 99.12.19:  Local 2D arrays trimmed!
!   WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!   WS  00.07.13:  Changed PILGRIM API
!   AAM 00.08.29:  Added kfirst, klast
!   AAM 01.06.27:  Added y communicators
!   SJL 30.07.01:  MPI optimization/simplification
!   WS  02.04.24:  New mod_comm interfaces
!   WS  02.07.04:  Fixed 2D decomposition bug dest/src for mp_send3d
!   WS  03.11.19:  Merged in CAM changes by Mirin
!   WS  03.12.03:  Added GRID as argument, dynamics_vars removed
!   WS  04.08.25:  Simplification of interface with GRID
!   WS  04.10.07:  Removed dependency on spmd_dyn; info now in GRID
!   WS  05.04.04:  Transitioned to type T_TRACERS (supports r4 and r8)
!   WS  05.04.09:  Each tracer now ghosted individually (save buffers)
!   WS  05.04.12:  Full support for either r4 or r8 tracers
!   WS  05.05.25:  Merged CAM and GEOS5, e.g. nsplt(k) opt. in CAM
!   PW  05.10.12:  Changes for Cray X1(E), alternative implementation
!                  of double buffering logic
!   WS  06.09.08:  Magic numbers are now F90 parameters
!
!EOP
!---------------------------------------------------------------------
!BOC

   real(r8), parameter ::  D1EM10                  =  1.0e-10_r8
   real(r8), parameter ::  D1_0                    =  1.0_r8
   real(r8), parameter ::  D0_0                    =  0.0_r8

! Local variables:
! 2d arrays
   real(r8)  a2(grid%im,grid%jfirst:grid%jlast)
   real(r8)  fx(grid%im,grid%jfirst:grid%jlast)
   real(r8)  fy(grid%im,grid%jfirst:grid%jlast+1)
   real(r8)  cymax(grid%kfirst:grid%klast)
! Temporary r8 array for Q
   real(r8)  ::  &
      q_r8(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d,grid%kfirst:grid%klast,1:2)

   real(r8) dp2(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
   logical ffsl(grid%jm,grid%kfirst:grid%klast)
   integer :: nsplt(grid%kfirst:grid%klast)

   integer :: im, jm, km                    ! Dimensions
   integer :: ng                            ! Max number of ghost latitudes
   integer :: jfirst, jlast, kfirst, klast  ! YZ decomposition limits
   integer :: cur, nxt                      ! current and next q_r8 buffer indices

   integer i, j, k
   integer it, iq, kq, max_nsplt
   integer :: k_courant, kend
   integer ktot
   integer js1gd, js2g0, js2gd, jn2g0,jn2gd,jn1g1,jn1gd
#if defined( SPMD )
   integer :: dest, src
#endif

   real(r8) cy_global
   real(r8) frac
   real(r8) cmax
   real(r8) sum1, sum2

   cur     = 1
   nxt     = 2

   im      = grid%im
   jm      = grid%jm
   km      = grid%km
   ng      = grid%ng_d

   jfirst  = grid%jfirst
   jlast   = grid%jlast
   kfirst  = grid%kfirst
   klast   = grid%klast

   ktot  = klast - kfirst + 1
   js2g0 = max(2,jfirst)
   jn2g0 = min(jm-1,jlast)
   jn1g1 = min(jm,jlast+1)
   js1gd = max(1,jfirst-ng)     ! NG latitudes on S (starting at 1)
   js2gd = max(2,jfirst-ng)     ! NG latitudes on S (starting at 2)
   jn2gd = min(jm-1,jlast+ng)   ! NG latitudes on S (ending at jm-1)
   jn1gd = min(jm,jlast+ng)     ! NG latitudes on N (ending at jm)

#if defined( SPMD )
      call FVstartclock(grid,'---TRAC2D_COMM')
      call mp_send4d_ns( grid%commyz, im, jm, km,                      &
                         1, jfirst, jlast, kfirst, klast,             &
                         ng, ng, cx )
! Send one latitude of both cy and mfy to the south
      dest = grid%iam-1
      src  = grid%iam+1
      if ( mod(grid%iam,grid%npr_y) == 0 ) dest = -1
      if ( mod(grid%iam+1,grid%npr_y) == 0 ) src = -1
      call mp_send3d_2( grid%commyz, dest, src, im, jm, km,            &
                        1, im, jfirst, jlast+1, kfirst, klast,        &
                        1, im, jfirst, jfirst, kfirst, klast, cy, mfy)
      call FVstopclock(grid,'---TRAC2D_COMM')
#endif

!$omp parallel do default(shared) private(i,j,k,cmax)
   do k=kfirst,klast
        cymax(k) = D0_0
       do j=js2g0,jlast
            cmax = D0_0
          do i=1,im
            cmax = max( abs(cy(i,j,k)), cmax)
          enddo
            cymax(k) = max(cymax(k), cmax*(D1_0 + grid%sine(j)**16) )
       enddo
   enddo

#if defined( SPMD )
   call FVstartclock(grid,'---TRAC2D_COMM')
   call mp_recv4d_ns( grid%commyz, im, jm, km,                         &
                      1, jfirst, jlast, kfirst, klast,                &
                      ng, ng, cx )
   call mp_recv3d_2( grid%commyz, src, im, jm, km,                     &
                     1, im, jfirst, jlast+1, kfirst, klast,           &
                     1, im, jlast+1, jlast+1, kfirst, klast, cy, mfy)

   call parcollective( grid%comm_y, MAXOP, ktot, cymax )
   call FVstopclock(grid,'---TRAC2D_COMM')
#endif

!---------------------------------------------------------------------
! Determine the required value of nsplt for each level
!---------------------------------------------------------------------
    nsplt(:)  = int( D1_0 + cymax(:) )
    max_nsplt = maxval( nsplt(:) )
#if defined( SPMD )
    call FVstartclock(grid,'---TRAC2D_COMM')
    call parcollective( grid%comm_z, MAXOP, max_nsplt )  ! Find global max
    call FVstopclock(grid,'---TRAC2D_COMM')
#endif
#ifndef WACCM_MOZART
    nsplt(:)  = max_nsplt
#endif
    do k_courant = klast,kfirst,-1
       if( nsplt(k_courant) > 1 ) then
          exit
       end if
    end do
    k_courant = max( kfirst,k_courant )
!!!    if (max_nsplt /= 1) write(iulog,*) 'trac2d: max_nsplt,k_courant = ', max_nsplt,k_courant
!!!    write(iulog,*) "max_nsplt", max_nsplt, "k_cour", k_courant, "nsplt", nsplt(:)

!$omp  parallel do default(shared) private(i,j,k,frac) schedule(dynamic,1)

#if !defined(USE_OMP)
!CSD$ PARALLEL DO PRIVATE (I, J, K, FRAC)
#endif
 do 4000 k=kfirst,klast

    if( nsplt(k) .ne. 1 ) then
       frac = D1_0 / nsplt(k)
       do j=js2gd,jn2gd                  
          do i=1,im
            cx(i,j,k) =  cx(i,j,k) * frac      ! cx ghosted on N*ng S*ng
          enddo
       enddo

       do j=js2g0,jn2g0
          do i=1,im
            mfx(i,j,k) = mfx(i,j,k) * frac
          enddo
       enddo

       do j=js2g0,jn1g1                     
          do i=1,im
             cy(i,j,k) =  cy(i,j,k) * frac    ! cy ghosted on N
            mfy(i,j,k) = mfy(i,j,k) * frac    ! mfy ghosted on N
          enddo
       enddo
    endif

       do j=js2g0,jn2g0
          do i=1,im
             if(cy(i,j,k)*cy(i,j+1,k) > D0_0) then
                if( cy(i,j,k) > D0_0) then
                   va(i,j,k) = cy(i,j,k)
                else
                   va(i,j,k) = cy(i,j+1,k)      ! cy ghosted on N
                endif
             else
                va(i,j,k) = D0_0
             endif
          enddo
       enddo

! Check if FFSL extension is needed.

       do j=js2gd,jn2gd             ! flux needed on N*ng S*ng
          ffsl(j,k) = .false.
          do i=1,im
            if( abs(cx(i,j,k)) > D1_0 ) then  ! cx ghosted on N*ng S*ng
              ffsl(j,k) = .true.
              exit
            endif
          enddo
       enddo

! Scale E-W mass fluxes by CX if FFSL
       do j=js2g0,jn2g0
          if( ffsl(j,k) ) then
            do i=1,im
              flx(i,j,k) = mfx(i,j,k) / sign( max(abs(cx(i,j,k)), D1EM10), &
                                        cx(i,j,k) )
            enddo
          else
            do i=1,im
              flx(i,j,k) = mfx(i,j,k)
            enddo
          endif
       enddo
4000  continue
#if !defined(USE_OMP)
!CSD$ END PARALLEL DO
#endif

 call FVbarrierclock(grid,'sync_trac2d_tracer',grid%commyz)
 call FVstartclock(grid,'---TRAC2D_TRACER')

 do 6000 it=1, max_nsplt
    if ( it == 1 ) then
       kend = klast       !  The entire vertical slab needs to be considered
    else
       kend = k_courant   !  Only the subset including courant # > 1 considered
    endif
! WS 05.04.06:  send only the first tracer the rest at end of do iq loop
!               NOTE: there is per definition at least one tracer
    q_r8(1:im,jfirst:jlast,kfirst:kend,1) = &
               tracer(1:im,jfirst:jlast,kfirst:kend,nlo)
#if defined( SPMD )
    call FVstartclock(grid,'---TRAC2D_TRACER_COMM')
    call mp_send4d_ns( grid%commyz, im, jm, km,                        &
                       1, jfirst, jlast, kfirst, kend,                &
                       ng, ng, q_r8(1,jfirst-ng,kfirst,1) )
    call FVstopclock(grid,'---TRAC2D_TRACER_COMM')
#endif

!$omp parallel do default(shared) private(i,j,k,sum1,sum2)

  do 3000 k=kfirst,kend
     if (it <= nsplt(k)) then
        do j=js2g0,jn2g0
           do i=1,im-1
              dp2(i,j,k) =  dp1(i,j,k) + mfx(i,j,k) - mfx(i+1,j,k) +  &
                        (mfy(i,j,k) - mfy(i,j+1,k)) * grid%acosp(j)
           enddo
           dp2(im,j,k) = dp1(im,j,k) + mfx(im,j,k) - mfx(1,j,k) +  &
                         (mfy(im,j,k) - mfy(im,j+1,k)) * grid%acosp(j)
        enddo

        if ( jfirst == 1  ) then
           sum1 = D0_0
           do i=1,im
              sum1 = sum1 + mfy(i,2,k)
           end do

           sum1 = - sum1 * grid%rcap
           do i=1,im
              dp2(i,1,k) = dp1(i,1,k) +  sum1
           enddo
        endif

        if ( jlast == jm ) then
           sum2 = D0_0
           do i=1,im
              sum2 = sum2 + mfy(i,jm,k)
           end do

              sum2 = sum2 * grid%rcap
           do i=1,im
              dp2(i,jm,k) = dp1(i,jm,k) +  sum2
           enddo
        endif
     endif
3000  continue

      do iq = nlo, nhi
#if defined( SPMD )
         call FVstartclock(grid,'---TRAC2D_TRACER_COMM')
!
! The buffer indices are exchanged, so that cur points to the current buffer,
! while nxt points to the one which will be used next.
!
        if ( mod(iq-nlo+1,2) == 0 ) then
          cur = 2
          nxt = 1
        else
          cur = 1
          nxt = 2
        endif
        call mp_recv4d_ns( grid%commyz, im, jm, km,                    &
                           1, jfirst, jlast, kfirst, kend,            &
                           ng, ng, q_r8(1,jfirst-ng,kfirst,cur) )

!
! Pre-send the next tracer
!
        if ( iq < nhi ) then
          q_r8(1:im,jfirst:jlast,kfirst:kend,nxt) = &
             tracer(1:im,jfirst:jlast,kfirst:kend,iq+1)
          call mp_send4d_ns( grid%commyz, im, jm, km,                  &
                             1, jfirst, jlast, kfirst, kend,          &
                             ng, ng, q_r8(1,jfirst-ng,kfirst,nxt) )
        endif
        call FVstopclock(grid,'---TRAC2D_TRACER_COMM')
#else
!
! No message passing -- simply copy the tracer into q_r8
!
        q_r8(1:im,jfirst:jlast,kfirst:kend,cur) = &
              tracer(1:im,jfirst:jlast,kfirst:kend,iq)
#endif

#if !defined(INNER_OMP)
!$omp parallel do default(shared)    &
!$omp private(i, j, k, kq, fx, fy, a2)
#endif
#if (!defined USE_OMP) 
!CSD$ PARALLEL DO PRIVATE (I, J, K, KQ, FX, FY, A2)
#endif
        do 5000 k=kfirst,kend
           if ( it <= nsplt(k) ) then
              call tp2c(a2, va(1,jfirst,k), q_r8(1:,jfirst-ng:,k,cur),  &
                        cx(1,jfirst-ng,k) , cy(1,jfirst,k),           &
                        im, jm, iord, jord, ng,                       &
                        fx, fy, ffsl(1,k), grid%rcap, grid%acosp,     &
                        flx(1,jfirst,k), mfy(1,jfirst,k),             &
                        grid%cosp, 1, jfirst, jlast )

              do j=jfirst,jlast
                 do i=1,im
                    q_r8(i,j,k,cur) = q_r8(i,j,k,cur)*dp1(i,j,k) + a2(i,j)
                 enddo
              enddo

              if (fill) call fillxy (q_r8(1:,jfirst:,k,cur), im, jm, jfirst, &
                                     jlast, grid%acap, grid%cosp, grid%acosp)

              do j=jfirst,jlast
                 do i=1,im
                    tracer(i,j,k,iq) = q_r8(i,j,k,cur) / dp2(i,j,k)
                 enddo
              enddo
           endif
5000    continue
#if (!defined USE_OMP) 
!CSD$ END PARALLEL DO
#endif

      enddo  ! End of do iq=nlo, nhi

!$omp parallel do private(i, j, k) schedule( dynamic,1 )
      do k=kfirst,kend
         if ( it < nsplt(k) ) then
            do j=jfirst,jlast
               do i=1,im
                  dp1(i,j,k) = dp2(i,j,k)
               enddo
            enddo
         endif
      enddo

6000  continue
 call FVstopclock(grid,'---TRAC2D_TRACER')

      return
!EOC
 end subroutine trac2d
!-----------------------------------------------------------------------


