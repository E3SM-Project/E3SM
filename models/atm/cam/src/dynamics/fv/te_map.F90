!-----------------------------------------------------------------------
!BOP
! !ROUTINE: te_map --- Map vertical Lagrangian coordinates to normal grid
!
! !INTERFACE:

   subroutine te_map(grid,    consv,  convt,   ps,      omga,            &
                     pe,      delp,    pkz,    pk,      mdt,             &
                     nx,      u,       v,      pt,      tracer,          &
                     hs,      cp,      akap,   kord,    peln,            &
                     te0,     te,      dz,     mfx,     mfy,             &
                     te_method )
!
! !USES:

   use shr_kind_mod,  only : r8 => shr_kind_r8
   use dynamics_vars, only : T_FVDYCORE_GRID
   use mapz_module,   only : map1_cubic_te, map1_ppm, mapn_ppm_tracer
   use cam_logfile,   only : iulog
   use abortutils,    only : endrun

#if defined( SPMD )
   use mod_comm,      only : mp_send3d, mp_recv3d
#endif
   use phys_control,  only: waccmx_is !WACCM-X runtime switch

   implicit none

#if defined( SPMD )
#define CPP_PRT_PREFIX  if(grid%iam==0)
#else
#define CPP_PRT_PREFIX
#endif

! !INPUT PARAMETERS:
   type (T_FVDYCORE_GRID), intent(inout) :: grid    ! grid for XY decomp
   logical consv                 ! flag to force TE conservation
   logical convt                 ! flag to control pt output (see below)
   integer mdt                   ! mapping time step (same as phys)
   integer nx                    ! number of SMP "decomposition" in x
   real(r8) hs(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy) ! surface geopotential
   real(r8) cp
   real(r8) te0
   integer, intent(in)          ::  te_method   ! Method for vertical total energy remapping
                                                !     0          :  piecewise-parabolic
                                                !     1          :  cubic interpolation
 
! !INPUT/OUTPUT PARAMETERS:
   real(r8) pk(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km+1) ! pe to the kappa
   real(r8) u(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)    ! u-wind (m/s)
   real(r8) v(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)    ! v-wind (m/s)
! tracers including specific humidity
!!!   real(r8) q(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km,grid%ntotq)
   real(r8), intent(inout) ::   &
       tracer(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km,grid%ntotq) ! Tracer array
   real(r8) pe(grid%ifirstxy:grid%ilastxy,grid%km+1,grid%jfirstxy:grid%jlastxy) ! pressure at layer edges
   real(r8) ps(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)      ! surface pressure
   real(r8) pt(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)   ! virtual potential temperature as input
                                     ! Output: virtual temperature if convt is true
                                     ! false: output is (virtual) potential temperature 
   real(r8)  te(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)  ! Work array (cache performance)
   real(r8)  dz(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)     ! Work array (cache performance)
   real(r8)  mfx(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)   
   real(r8)  mfy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)  
 
! !OUTPUT PARAMETERS:
   real(r8) delp(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)    ! pressure thickness
   real(r8) omga(grid%ifirstxy:grid%ilastxy,grid%km,grid%jfirstxy:grid%jlastxy)    ! vertical press. velocity (pascal/sec)
   real(r8) peln(grid%ifirstxy:grid%ilastxy,grid%km+1,grid%jfirstxy:grid%jlastxy)  ! log(pe)
   real(r8) pkz(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)     ! layer-mean pk for converting t to pt

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! WS 99.05.19 : Replaced IMR, JMR, JNP and NL with IM, JM-1, JM and KM
! WS 99.05.25 : Revised conversions with IMR and JMR; removed fvcore.h
! WS 99.07.29 : Reduced computation region to grid%jfirstxy:jlast
! WS 99.07.30 : Tested concept with message concatenation of te_map calls
! WS 99.10.01 : Documentation; indentation; cleaning
! SJL 99.12.31: SMP "decomposition" in-E-W direction
! WS 00.05.14 : Renamed ghost indices as per Kevin's definitions
! WS 00.07.13 : Changed PILGRIM API
! AM 00.08.29 : Variables in this routine will ultimately be decomposed in (i,j).
! AM 01.06.13 : 2-D decomposition; reordering summation causes roundoff difference.
! WS 01.06.10 : Removed "if(first)" section in favor of a variable module
! AM 01.06.27 : Merged yz decomposition technology into ccm code.
! WS 02.01.14 : Upgraded to mod_comm
! WS 02.04.22 : Use mapz_module from FVGCM
! WS 02.04.25 : New mod_comm interfaces
! WS 03.08.12 : Introduced unorth
! WS 03.11.19 : Merged in CAM changes by Mirin
! WS 03.12.03 : Added GRID as argument, dynamics_vars removed
! WS 04.08.25 : Simplified interface by using GRID
! WS 04.10.07 : Removed dependency on spmd_dyn; info now in GRID 
! WS 05.03.25 : Changed tracer to type T_TRACERS
! WS 05.04.12 : Call mapn_ppm_tracer instead of mapn_ppm
! AT 05.05.11 : Merged with the version Cerebus (unique pole issues)
! WS 05.05.18 : Merged CAM and GEOS5 versions (mostly GEOS5)
! LT 05.11.14 : Call map1_cubic_te for Cubic Interpolation of Total Energy
! WP 06.01.18 : Added calls to map1_ppm for horizontal mass fluxes
! LT 06.02.08 : Implement code for partial remapping option
! WS 06.11.29 : Merge CAM/GEOS5; magic numbers isolated; te_method
! CC 07.01.29 : Additions for proper calculation of OMGA
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

! Magic numbers used in this module
      real(r8), parameter ::  D0_0                    =  0.0_r8
      real(r8), parameter ::  D0_25                   =  0.25_r8
      real(r8), parameter ::  D0_5                    =  0.5_r8
      real(r8), parameter ::  D1_0                    =  1.0_r8
      real(r8), parameter ::  D2_0                    =  2.0_r8
      real(r8), parameter ::  D10_0                   = 10.0_r8
      real(r8), parameter ::  D1E4                    =  1.0e4_r8

      integer :: im, jm, km            ! x, y, z dimensions
      integer :: nq                    ! number of tracers to be advected
      integer :: ntotq                 ! Total number of tracers
      integer :: ifirst, ilast         ! starting & ending longitude index
      integer :: jfirst, jlast         ! starting & ending latitude index
      integer :: myidxy_y, iam
      integer :: nprxy_x, nprxy_y

! Local variables for Partial Remapping
! -------------------------------------
      real(r8) ::    pref(grid%km+1)
      real(r8) ::     fac(grid%km+1)
      real(r8) ::      zz(grid%km+1)
      real(r8) ::      z1,z2

      real(r8), parameter :: alf   =    0.042_r8
      real(r8), parameter :: pa    =      1.0_r8
      real(r8), parameter :: pb    =    500.0_r8
      real(r8), parameter :: psurf = 100001.0_r8
      real(r8), parameter :: bet   = D2_0*alf/(D1_0+alf)

! Local arrays:
! -------------
      real(r8) rmin(nx*grid%jm), rmax(nx*grid%jm)
      real(r8) tte(grid%jm)
! x-y
      real(r8)  u2(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy+1)
      real(r8)  v2(grid%ifirstxy:grid%ilastxy+1,grid%jfirstxy:grid%jlastxy)
      real(r8)  t2(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)
      real(r8)  veast(grid%jfirstxy:grid%jlastxy,grid%km)
! y-z
      real(r8)  pe0(grid%ifirstxy:grid%ilastxy,grid%km+1)
      real(r8)  pe1(grid%ifirstxy:grid%ilastxy,grid%km+1)
      real(r8)  pe2(grid%ifirstxy:grid%ilastxy,grid%km+1)
      real(r8)  pe3(grid%ifirstxy:grid%ilastxy,grid%km+1)
      real(r8) phis(grid%ifirstxy:grid%ilastxy,grid%km+1)
      real(r8) u2_sp(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) v2_sp(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) t2_sp(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) u2_np(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) v2_np(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) t2_np(grid%ifirstxy:grid%ilastxy,grid%km)

! x
      real(r8)     gz(grid%ifirstxy:grid%ilastxy)
      real(r8)  ratio(grid%ifirstxy:grid%ilastxy)
      real(r8)    bte(grid%ifirstxy:grid%ilastxy)
! z
      real(r8) pe1w(grid%km+1)
      real(r8) pe2w(grid%km+1)

      integer i1w, nxu
      integer i, j, k, js2g0, jn2g0, jn1g1
      integer kord
      integer krd

      real(r8) akap, dak, bkh, qmax, qmin
      real(r8) te_sp(grid%km), te_np(grid%km)
      real(r8) xysum(grid%jfirstxy:grid%jlastxy,2)
      real(r8) tmpik(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) tmpij(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,2)
      real(r8) omga_ik(grid%ifirstxy:grid%ilastxy,grid%km)    ! vertical press. velocity (tmp 2-d array)
      real(r8) dtmp
      real(r8) sum
      real(r8) rdt5
      real(r8) rg
      real(r8) te1
      real(r8) dlnp
      real(r8) tvm

      integer ixj, jp, it, i1, i2

#if defined( SPMD )
      integer :: dest, src
      real(r8)  unorth(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8)  pewest(grid%km+1,grid%jfirstxy:grid%jlastxy)
      real(r8), allocatable :: pesouth(:,:)
#endif

      integer comm_use, npry_use, itot

      logical diag

      data diag    /.false./

      z1 = log(pa/psurf)
      z2 = log(pb/psurf)

      im = grid%im
      jm = grid%jm
      km = grid%km
      nq = grid%nq

      ifirst = grid%ifirstxy
      ilast  = grid%ilastxy
      jfirst = grid%jfirstxy
      jlast  = grid%jlastxy

      iam      = grid%iam
      myidxy_y = grid%myidxy_y
      nprxy_x  = grid%nprxy_x
      nprxy_y  = grid%nprxy_y

! Intialize PREF and FAC for Partial Remapping (above 100-mb)
! -----------------------------------------------------------
      do k=1,km+1
         pref(k) = grid%ak(k) + grid%bk(k)*psurf
      enddo
      zz  = log( pref/psurf )
      zz  = D10_0*(zz-z2)/z1
      zz  = (D1_0-bet)*tanh(zz)
      do k=1,km+1
!     if( pref(k).le.D1E4 ) then 
!          fac(k) = (D1_0-zz(k))/(D2_0-bet)
!     else
           fac(k) = D1_0
!     endif
      enddo

! WS 99.07.27 :  Set loop limits appropriately
! --------------------------------------------
      js2g0  = max(2,jfirst)
      jn1g1  = min(jm,jlast+1)
      jn2g0  = min(jm-1,jlast)
      do j=jfirst,jlast
         xysum(j,1) = D0_0
         xysum(j,2) = D0_0
      enddo
      do j=jfirst,jlast
         do i=ifirst,ilast
            tmpij(i,j,1) = D0_0
            tmpij(i,j,2) = D0_0
         enddo
      enddo
      do k=1,km
         do i=ifirst,ilast
            tmpik(i,k) = D0_0
         enddo
      enddo

      itot = ilast-ifirst+1
      nxu = 1
      if (itot == im) nxu = nx

#if defined( SPMD )
      comm_use = grid%commxy_y
      npry_use = nprxy_y

      call mp_send3d( grid%commxy, iam-nprxy_x, iam+nprxy_x, im, jm, km,  &
                      ifirst, ilast, jfirst, jlast, 1, km,               &
                      ifirst, ilast, jfirst, jfirst, 1, km, u )
! Nontrivial x decomposition
      if (itot /= im) then
        dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
        src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
        call mp_send3d( grid%commxy, dest, src, im, jm, km,               &
                        ifirst, ilast, jfirst, jlast, 1,km,              &
                        ifirst, ifirst, jfirst, jlast, 1, km, v )
      endif
#endif
      call pkez(nxu, im, km, jfirst, jlast, 1, km, ifirst, ilast,        &
                pe, pk, akap, grid%ks, peln, pkz, .false.)
 
! Single subdomain case (periodic)
      do k=1,km
         do j=jfirst,jlast
            veast(j,k) = v(ifirst,j,k)
         enddo
      enddo
#if defined( SPMD ) 
      call mp_recv3d( grid%commxy, iam+nprxy_x, im, jm, km,               &
                      ifirst, ilast, jlast+1, jlast+1, 1, km,            &
                      ifirst, ilast, jlast+1, jlast+1, 1, km, unorth )
! Nontrivial x decomposition
      if (itot /= im) then
        call mp_recv3d( grid%commxy, src, im, jm, km,                     &
                        ilast+1, ilast+1, jfirst, jlast, 1, km,          &
                        ilast+1, ilast+1, jfirst, jlast, 1, km, veast )
        dest = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
        src  = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
        call mp_send3d( grid%commxy, dest, src, im, km+1, jm,             &
                        ifirst, ilast, 1, km+1, jfirst, jlast,           &
                        ilast, ilast, 1, km+1, jfirst, jlast, pe )
      endif
      call mp_send3d( grid%commxy, iam+nprxy_x, iam-nprxy_x, im, km+1,jm, &
                      ifirst, ilast, 1, km+1, jfirst, jlast,             &
                      ifirst, ilast, 1, km+1, jlast, jlast, pe )
#endif

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(i,j, k, u2, v2, t2)

! Compute cp*T + KE

      do 1000 k=1,km

        do j=js2g0,jlast
           do i=ifirst,ilast
              u2(i,j) = u(i,j,k)**2
           enddo
        enddo
#if defined( SPMD )
        if ( jlast < jm ) then
           do i=ifirst,ilast
              u2(i,jlast+1) = unorth(i,k)**2  ! fill ghost zone
           enddo
        endif
#endif

        do j=js2g0,jn2g0
          do i=ifirst,ilast
             v2(i,j) = v(i,j,k)**2
          enddo
          v2(ilast+1,j) = veast(j,k)**2
        enddo

        do j=jfirst,jlast
           do i=ifirst,ilast
              t2(i,j) = cp*pt(i,j,k)
           enddo
        enddo

        do j=js2g0,jn2g0
          do i=ifirst,ilast
            te(i,j,k) = D0_25 * ( u2(i,j) + u2(i,j+1) +        &
                                  v2(i,j) + v2(i+1,j) ) +      &
                        t2(i,j)*pkz(i,j,k)
          enddo
        enddo

! WS 99.07.29 : Restructuring creates small round-off.  Not clear why...

! Do collective Mpisum (in i) for te_sp and te_np below (AAM)
!
        if ( jfirst == 1 ) then
! South pole
          do i=ifirst,ilast
             u2_sp(i,k) = u2(i,2)
             v2_sp(i,k) = v2(i,2)
             t2_sp(i,k) = t2(i,1)
          enddo
        endif

        if ( jlast == jm ) then
! North pole
          do i=ifirst,ilast
             u2_np(i,k) = u2(i,jm)
             v2_np(i,k) = v2(i,jm-1)
             t2_np(i,k) = t2(i,jm)
          enddo
        endif

! Compute dz; geo-potential increments
        do j=jfirst,jlast
           do i=ifirst,ilast
              dz(i,j,k) = t2(i,j)*(pk(i,j,k+1)-pk(i,j,k))
           enddo
        enddo
1000  continue

#if defined( SPMD )
      allocate( pesouth(ifirst:ilast,km+1) )
      if (itot /= im) then
        call mp_recv3d( grid%commxy, src, im, km+1, jm,                   &
                        ifirst-1, ifirst-1, 1, km+1, jfirst, jlast,      &
                        ifirst-1, ifirst-1, 1, km+1, jfirst, jlast, pewest )   
      endif
      call mp_recv3d( grid%commxy, iam-nprxy_x, im, km+1, jm,             &
                      ifirst, ilast, 1, km+1, jfirst-1, jfirst-1,        &
                      ifirst, ilast, 1, km+1, jfirst-1, jfirst-1, pesouth )
#endif

      if ( jfirst == 1 ) then

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(i, k)

         do k = 1, km
            do i=ifirst,ilast
              tmpik(i,k) = D0_5*( u2_sp(i,k) + v2_sp(i,k) ) + t2_sp(i,k)*pkz(i,1,k)
            enddo
         enddo

         call par_xsum( grid, tmpik, km, te_sp)

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(i, k)

         do k = 1, km
            te_sp(k) = te_sp(k)/real(im,r8)
            do i=ifirst,ilast
              te(i,  1,k) = te_sp(k)
            enddo
         enddo
      endif

      if ( jlast == jm ) then

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            do i=ifirst,ilast
              tmpik(i,k) = D0_5*( u2_np(i,k) + v2_np(i,k) ) + t2_np(i,k)*pkz(i,jm,k)
            enddo
         enddo

         call par_xsum( grid, tmpik, km, te_np)

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            te_np(k) = te_np(k)/real(im,r8)
            do i=ifirst,ilast
              te(i,jm,k) = te_np(k)
            enddo
         enddo
      endif

      it = itot / nxu
      jp = nxu * ( jlast - jfirst + 1 )

!$omp  parallel do           &
!$omp  default(shared)       &
!$omp  private(i,j,k,i1w,pe0,pe1,pe2,pe3,ratio)   &
!$omp  private(dak,bkh,rdt5,phis,krd, ixj,i1,i2) &
!$omp  private(pe1w, pe2w, omga_ik )

!     do 2000 j=jfirst,jlast
      do 2000 ixj=1,jp

        j  = jfirst + (ixj-1) / nxu
        i1 = ifirst + it * mod(ixj-1, nxu)
        i2 = i1 + it - 1

! Copy data to local 2D arrays.
        i1w = i1-1
        if (i1 == 1) i1w = im
        do k=1,km+1
           do i=i1,i2
              pe1(i,k) = pe(i,k,j)
              if (k>1) then
                if (pe1(i,k)-pe1(i,k-1)<1.0e-7_r8) then
                  write(iulog,*) "Lagrangian levels are crossing"
                  write(iulog,*) "Run will ABORT!"
                  write(iulog,*) "Suggest to increase NSPLTVRM"
                  call endrun
                endif
              endif
           enddo
           if( itot == im ) then
               pe1w(k) = pe(i1w,k,j)
#if defined( SPMD )
           else
               pe1w(k) = pewest(k,j)
#endif
           endif
        enddo

        do k=1,grid%ks+1
           do i=i1,i2
              pe0(i,k) = grid%ak(k)
              pe2(i,k) = grid%ak(k)
              pe3(i,k) = grid%ak(k)
            enddo
        enddo

        do k=grid%ks+2,km
           do i=i1,i2
              pe0(i,k) = grid%ak(k) + grid%bk(k)* ps(i,j)     ! Remapped PLE based on Old     PS
              pe2(i,k) = grid%ak(k) + grid%bk(k)*pe1(i,km+1)  ! Remapped PLE based on Updated PS
           enddo
        enddo

        do i=i1,i2
           pe0(i,km+1) =  ps(i,j)
           pe2(i,km+1) = pe1(i,km+1)
        enddo

! Ghosting for v mapping
        do k=grid%ks+2,km
           pe2w(k) = grid%ak(k) + grid%bk(k)*pe1w(km+1)
        enddo
        pe2w(km+1) = pe1w(km+1)

! update ps
! ---------
        do i=i1,i2
          ps(i,j)   = pe1(i,km+1)
        enddo

! Compute GeoPotential Heights and add to Total Energy
! Note:   GeoPotential Height(L) => d(Pe*PHI)/dPe
! ----------------------------------------------------

        do i=i1,i2
           phis(i,km+1) = hs(i,j)      
        enddo

        do k=km,1,-1
          do i=i1,i2
             phis(i,k) = phis(i,k+1) + dz(i,j,k)   
          enddo
        enddo

        do k=1,km+1
          do i=i1,i2
             phis(i,k) = phis(i,k) * pe1(i,k)
          enddo
        enddo

        do k=1,km
          do i=i1,i2
            te(i,j,k) =  te(i,j,k) + (phis(i,k+1) - phis(i,k)) /   &
                         (pe1(i,k+1) - pe1(i,k) )
          enddo
        enddo

! #######################################################################
! #                           ReMap Total Energy
! #######################################################################

! Compute Target Pressures for Partial Remapping
! ----------------------------------------------
        do k=1,km+1
          do i=i1,i2
            pe2(i,k) = fac(k)*pe2(i,k) + (D1_0-fac(k))*pe1(i,k)
          enddo
        enddo

! Update Delta-Pressure (from final remapped pressures)
! -----------------------------------------------------
        do k=1,km
          do i=i1,i2
             delp(i,j,k) = pe2(i,k+1) - pe2(i,k)
          enddo
        enddo

! Compute omga (dp/dt)
! --------------------
        rdt5 = D0_5 / real(mdt,r8)
        do k=2,km+1
           do i=i1,i2
              pe0(i,k) = pe1(i,k) - pe0(i,k)  ! Delta-P:  PLE(After CD_Core) minus PLE(Remapped based on old PS)
           enddo
        enddo

! ReMap Total Energy
! ------------------
        select case ( te_method )
        case (0)
          call map1_ppm ( km,   pe1,   te,                      &
                          km,   pe2,   te,  0,  0,              &
                          itot, i1-ifirst+1, i2-ifirst+1,       &
                          j, jfirst, jlast, 1, kord)
        case (1)
          call map1_cubic_te ( km,   pe1,   te,                 &
                               km,   pe2,   te,  0,  0,         &
                               itot, i1-ifirst+1, i2-ifirst+1,  &
                               j, jfirst, jlast, 1, kord)
        case default
          call map1_ppm ( km,   pe1,   te,                      &
                          km,   pe2,   te,  0,  0,              &
                          itot, i1-ifirst+1, i2-ifirst+1,       &
                          j, jfirst, jlast, 1, kord)
        endselect

! C.-C. Chen
! Map omga
        do k=1,km
           do i=i1,i2
              omga_ik(i,k) = omga(i,k,j)

              if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
                if (k == 1) omga_ik(i,k) = D0_0  
              endif

           end do
        end do
        call map1_ppm ( km,   pe1,   omga_ik,                   &
                        km,   pe2,   omga_ik,  0,  0,           &
                        itot, i1-ifirst+1, i2-ifirst+1,         &
                        1,      1,     1, 1, kord)
        do k=1,km
           do i=i1,i2
              omga(i,k,j) = omga_ik(i,k)
           end do
        end do

! #######################################################################
! #                           ReMap Constituents
! #######################################################################

       if( nq /= 0 ) then
          if(kord == 8) then
             krd = 8
          else
             krd = 7
          endif

          call mapn_ppm_tracer ( km,   pe1,   tracer, nq,                &
                                 km,   pe2,   i1, i2,                    &
                                 j, ifirst, ilast, jfirst, jlast, 0, krd)
       endif

! #######################################################################
! #                             ReMap U-Wind
! #######################################################################

        if(j /= 1) then

! WS 99.07.29 : protect j==jfirst case
          if (j > jfirst) then
            do k=2,km+1
              do i=i1,i2
                pe0(i,k) = D0_5*(pe1(i,k)+pe(i,k,j-1))
              enddo
            enddo

            do k=grid%ks+2,km+1
              bkh = D0_5*grid%bk(k)
              do i=i1,i2
                pe3(i,k) = grid%ak(k) + bkh*(pe1(i,km+1)+pe(i,km+1,j-1))
              enddo
            enddo

#if defined( SPMD )
          else
!  WS 99.10.01 : Read in pe(:,:,jfirst-1) from the pesouth buffer
            do k=2,km+1
              do i=i1,i2
                pe0(i,k) = D0_5*(pe1(i,k)+pesouth(i,k))
              enddo
            enddo

            do k=grid%ks+2,km+1
              bkh = D0_5*grid%bk(k)
              do i=i1,i2
                pe3(i,k) = grid%ak(k) + bkh*(pe1(i,km+1)+pesouth(i,km+1))
              enddo
            enddo
#endif
          endif


! Compute Target Pressures for Partial Remapping
! ----------------------------------------------
          do k=1,km+1
            do i=i1,i2
              pe3(i,k) = fac(k)*pe3(i,k) + (D1_0-fac(k))*pe0(i,k)
            enddo
          enddo

! ReMap U-Wind (D-Grid Location)
! ------------------------------
          call map1_ppm ( km,   pe0,    u,    km,   pe3,    u,            &
                          0,    0,   itot, i1-ifirst+1, i2-ifirst+1,      &
                          j,    jfirst, jlast,  -1,    kord)


! ReMap Y-Mass Flux (C-Grid Location)
! -----------------------------------
          do k=1,km
             do i=i1,i2
                mfy(i,j,k) = mfy(i,j,k)/(pe0(i,k+1)-pe0(i,k))
             enddo
          enddo
          call map1_ppm ( km,   pe0,    mfy,    km,   pe3,    mfy,        &
                          0,    0,   itot, i1-ifirst+1, i2-ifirst+1,      &
                          j,    jfirst, jlast,  -1,    kord)
          do k=1,km
             do i=i1,i2
                mfy(i,j,k) = mfy(i,j,k)*(pe3(i,k+1)-pe3(i,k))
             enddo
          enddo
        endif

! #######################################################################
! #                             ReMap V-Wind
! #######################################################################

        if(j /= 1 .and. j /= jm) then
          do k=2,km+1
! pe1(i1-1,1:km+1) must be ghosted
            pe0(i1,k) = D0_5*(pe1(i1,k)+pe1w(k))
            do i=i1+1,i2
               pe0(i ,k) = D0_5*(pe1(i,k)+pe1(i-1,k))
            enddo
          enddo

          do k=grid%ks+2,km+1
! pe2(i1-1,grid%ks+2:km+1) must be ghosted
            pe3(i1,k) = D0_5*(pe2(i1,k)+pe2w(k))
            do i=i1+1,i2
               pe3(i,k) = D0_5*(pe2(i,k)+pe2(i-1,k))
            enddo
          enddo

! Compute Target Pressures for Partial Remapping
! ----------------------------------------------
          do k=1,km+1
            do  i=i1,i2
              pe3(i,k) = fac(k)*pe3(i,k) + (D1_0-fac(k))*pe0(i,k)
            enddo
          enddo

! ReMap V-Wind (D-Grid Location)
! ------------------------------
          call map1_ppm ( km,   pe0,    v,    km,   pe3,    v,         &
                          0,    0,   itot, i1-ifirst+1, i2-ifirst+1,   &
                          j,    jfirst, jlast, -1, kord)

! ReMap X-Mass Flux (C-Grid Location)
! -----------------------------------
          do k=1,km
             do i=i1,i2
                mfx(i,j,k) = mfx(i,j,k)/(pe0(i,k+1)-pe0(i,k))
             enddo
          enddo
          call map1_ppm ( km,   pe0,   mfx,    km,   pe3,   mfx,       &
                          0,    0,   itot, i1-ifirst+1, i2-ifirst+1,   &
                          j, jfirst, jlast, -1, kord)
          do k=1,km
             do i=i1,i2
                mfx(i,j,k) = mfx(i,j,k)*(pe3(i,k+1)-pe3(i,k))
             enddo
          enddo
        endif

! Save new PE to temp storage peln
! --------------------------------
        do k=2,km
          do i=i1,i2
             peln(i,k,j) = pe2(i,k)
          enddo
        enddo

! Check deformation.
       if( diag ) then
          rmax(ixj) = D0_0
          rmin(ixj) = D1_0
          do k=1,km
             do i=i1,i2
              ratio(i) = (pe1(i,k+1)-pe1(i,k)) / (pe2(i,k+1)-pe2(i,k))
             enddo

             do i=i1,i2
              if(ratio(i) > rmax(ixj)) then
                 rmax(ixj) = ratio(i)
              elseif(ratio(i) < rmin(ixj)) then
                 rmin(ixj) = ratio(i)
              endif
            enddo
          enddo
       endif
2000  continue


#if defined( SPMD )
      deallocate( pesouth )

! Send u southward
      call mp_send3d( grid%commxy, iam-nprxy_x, iam+nprxy_x, im, jm, km,&
                      ifirst, ilast, jfirst, jlast, 1, km,             &
                      ifirst, ilast, jfirst, jfirst, 1, km, u )
      if (itot /= im) then
        dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
        src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
        call mp_send3d( grid%commxy, dest, src, im, jm, km,             &
                        ifirst, ilast, jfirst, jlast, 1, km,           &
                        ifirst, ifirst, jfirst, jlast, 1, km, v )
      endif
#endif

      if( diag ) then
        qmin = rmin(1)
        do ixj=2, jp
          if(rmin(ixj) < qmin) then
            qmin = rmin(ixj)
          endif
        enddo
        CPP_PRT_PREFIX write(iulog,*) 'rmin=', qmin

        qmax = rmax(1)
        do ixj=2, jp
          if(rmax(ixj) > qmax) then
            qmax = rmax(ixj)
          endif
        enddo
        CPP_PRT_PREFIX write(iulog,*) 'rmax=', qmax
      endif

! Recover Final Edge-Pressures and Compute Mid-Level PKZ
! ------------------------------------------------------

!$omp  parallel do          &
!$omp  default(shared)      &
!$omp  private(i,j,k)

      do j=jfirst,jlast
        do k=2,km
          do i=ifirst,ilast
            pe(i,k,j) = peln(i,k,j)
          enddo
        enddo
      enddo

      do k=1,km+1
        do j=jfirst,jlast
          do i=ifirst,ilast
            pk(i,j,k) = pe(i,k,j)**akap
          enddo
        enddo
      enddo
      call pkez(nxu, im, km, jfirst, jlast, 1, km, ifirst, ilast,  &
                pe, pk, akap, grid%ks, peln, pkz, .false.)

! Single x-subdomain case (periodic)
      do k = 1, km
        do j = jfirst, jlast
          veast(j,k) = v(ifirst,j,k)
        enddo
      enddo

#if defined( SPMD )
! Recv u from north
      call mp_recv3d( grid%commxy, iam+nprxy_x, im, jm, km,              &
                      ifirst, ilast, jlast+1, jlast+1, 1, km,           &
                      ifirst, ilast, jlast+1, jlast+1, 1, km, unorth )
      if (itot /= im) then
        call mp_recv3d( grid%commxy, src, im, jm, km,                    &
                        ilast+1, ilast+1, jfirst, jlast, 1, km,         &
                        ilast+1, ilast+1, jfirst, jlast, 1, km, veast )
      endif
#endif

! ((((((((((((((((( compute globally integrated TE )))))))))))))))))

      if( consv ) then

!$omp  parallel do         &
!$omp  default(shared)     &
!$omp  private(i,j,k)

        do k=1,km
          do j=jfirst,jlast
             do i=ifirst,ilast
                dz(i,j,k) = te(i,j,k) * delp(i,j,k)
             enddo
          enddo
        enddo

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(i,j,k,bte)

! Perform vertical integration

        do 4000 j=jfirst,jlast

          if ( j == 1 ) then
! SP
            tte(1) = D0_0
  
            do k=1,km
              tte(1) = tte(1) + dz(ifirst,1,k)
            enddo

          elseif ( j == jm) then
! NP
            tte(jm) = D0_0

            do k=1,km
              tte(jm) = tte(jm) + dz(ifirst,jm,k)
            enddo

          else
! Interior
            do i=ifirst,ilast
              bte(i) = D0_0
            enddo

            do k=1,km
              do i=ifirst,ilast
                bte(i) = bte(i) + dz(i,j,k)
              enddo
            enddo

            do i=ifirst,ilast
              tmpij(i,j,1) = bte(i)
            enddo

          endif
4000    continue

        call par_xsum( grid, tmpij, jlast-jfirst+1, xysum)

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(j)

        do j = max(jfirst,2), min(jlast,jm-1)
           tte(j) = xysum(j,1)*grid%cosp(j)
        enddo

        if ( jfirst == 1 ) tte(1)  = grid%acap * tte(1)
        if ( jlast == jm ) tte(jm) = grid%acap * tte(jm)

        te1 = D0_0
        call par_vecsum(jm, jfirst, jlast, tte, te1, comm_use, npry_use)

      endif   ! consv

      if( consv ) then

!$omp  parallel do       &
!$omp& default(shared)   &
!$omp& private(i,j)
 
       do j=js2g0, jn2g0
        do i=ifirst,ilast
          tmpij(i,j,1) = ps(i,j)
          tmpij(i,j,2) = peln(i,km+1,j) 
        enddo
       enddo

       call par_xsum( grid, tmpij, 2*(jlast-jfirst+1), xysum)

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(j)
 
       do j=js2g0, jn2g0
        tte(j) = cp*grid%cosp(j)*(xysum(j,1) - grid%ptop*real(im,r8) -           &
                 akap*grid%ptop*(xysum(j,2) - peln(ifirst,1,j)*real(im,r8)) )
! peln(i,1,j) should be independent of i (AAM)
       enddo

       if ( jfirst == 1 ) tte(1) = grid%acap*cp * (ps(ifirst,1) - grid%ptop -    &
                akap*grid%ptop*(peln(ifirst,km+1,1) - peln(ifirst,1,1) ) )
       if ( jlast == jm ) tte(jm)= grid%acap*cp * (ps(ifirst,jm) - grid%ptop -   &
                akap*grid%ptop*(peln(ifirst,km+1,jm) - peln(ifirst,1,jm) ) )
      endif ! consv

      if (consv) then

       sum=D0_0
       call par_vecsum(jm, jfirst, jlast, tte, sum, comm_use, npry_use)

       dtmp = (te0 - te1) / sum
       if( diag ) then
         CPP_PRT_PREFIX write(iulog,*) 'te=',te0, ' Energy deficit in T = ', dtmp
       endif

      endif              ! end consv check

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i,j,k, u2, v2)

! --------------------------------------------------------------------
! ---   Recover Tv from remapped Total Energy and its components   ---
! --------------------------------------------------------------------

      do 8000 k=1,km

! Intialize Kinetic Energy
! ------------------------
        do j=js2g0,jlast
          do i=ifirst,ilast
            u2(i,j) = u(i,j,k)**2
          enddo
        enddo
#if defined( SPMD )
        if ( jlast < jm ) then
           do i=ifirst,ilast
              u2(i,jlast+1) = unorth(i,k)**2  ! fill ghost zone
           enddo
        endif
#endif

        do j=js2g0,jn2g0
          do i=ifirst,ilast
            v2(i,j) = v(i,j,k)**2
          enddo
          v2(ilast+1,j) = veast(j,k)**2
        enddo

! Subtract Kinetic Energy from Total Energy (Leaving Internal + Potential)
! ------------------------------------------------------------------------
        do j=js2g0,jn2g0
          do i=ifirst,ilast
            te(i,j,k) = te(i,j,k) - D0_25 * ( u2(i,j) + u2(i,j+1)     &
                                             +v2(i,j) + v2(i+1,j) )
          enddo
        enddo

! South pole
! ----------
        if ( jfirst == 1 ) then
          do i=ifirst,ilast
            u2_sp(i,k) = u2(i,2)
            v2_sp(i,k) = v2(i,2)
          enddo
        endif

! North pole
! ----------
        if ( jlast == jm ) then
          do i=ifirst,ilast
            u2_np(i,k) = u2(i,jm)
            v2_np(i,k) = v2(i,jm-1)
          enddo
        endif

8000  continue

! South pole
! ----------
      if ( jfirst == 1 ) then

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            do i=ifirst,ilast
              tmpik(i,k) = D0_5*( u2_sp(i,k) + v2_sp(i,k) )
            enddo
         enddo

         call par_xsum( grid, tmpik, km, te_sp)

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            te_sp(k) = te_sp(k)/real(im,r8)
            do i=ifirst,ilast
              te(i,1,k) = te(i,1,k) - te_sp(k)
            enddo
         enddo
      endif

! North pole
! ----------
      if ( jlast == jm ) then

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            do i=ifirst,ilast
              tmpik(i,k) = D0_5*( u2_np(i,k) + v2_np(i,k) )
            enddo
         enddo

         call par_xsum( grid, tmpik, km, te_np)

!$omp  parallel do       &
!$omp  default(shared)   &
!$omp  private(i, k)

         do k = 1, km
            te_np(k) = te_np(k)/real(im,r8)
            do i=ifirst,ilast
              te(i,jm,k) = te(i,jm,k) - te_np(k)
            enddo
         enddo
      endif

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(ixj, i1, i2, i, j, k, rg, gz, tvm, dlnp)

! Recover Tv from remapped Total Energy and its components
! --------------------------------------------------------
      do 9000 ixj=1,jp

         j  = jfirst + (ixj-1) / nxu
         i1 = ifirst + it * mod(ixj-1, nxu)
         i2 = i1 + it - 1

         rg = akap * cp
         do i=i1,i2
            gz(i) = hs(i,j)           ! Initialize GeoPotential Heights
         enddo

        do k=km,1,-1
          do i=i1,i2
            dlnp  = rg*(peln(i,k+1,j) - peln(i,k,j))
            tvm   = delp(i,j,k)*(te(i,j,k) - gz(i)) /     &
                     ( cp*delp(i,j,k) - pe(i,k,j)*dlnp )
            gz(i) = gz(i) + dlnp*tvm  ! Update GeoPotential Heights
            pt(i,j,k) = tvm           ! pt is now (virtual) temperature
          enddo

          if( consv ) then
              do i=i1,i2
                 pt(i,j,k) = pt(i,j,k) + dtmp
              enddo
          endif

          if( .not. convt ) then
              do i=i1,i2
                 pt(i,j,k) = pt(i,j,k) / pkz(i,j,k)  ! Scaled Virtual Potential Temperature
              enddo
          endif
        enddo           ! end k-loop
9000  continue

      return
!EOC
      end subroutine te_map
!-----------------------------------------------------------------------
