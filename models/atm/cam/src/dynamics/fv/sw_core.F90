module sw_core
!BOP
!
! !MODULE: sw_core --- Utilities for solving the shallow-water equation
!
! !USES:
  use dynamics_vars, only: T_FVDYCORE_GRID
  use shr_kind_mod, only : r8 => shr_kind_r8

#ifdef NO_R16
   integer,parameter :: r16= selected_real_kind(12) ! 8 byte real
#else
   integer,parameter :: r16= selected_real_kind(24) ! 16 byte real
#endif

!
! !PUBLIC MEMBER FUNCTIONS:
      public d2a2c_winds, c_sw, d_sw
!
! !DESCRIPTION:
!
! This module contains vertical independent part of the Lagrangian
! dynamics; in simpler terms, it solves the 2D shallow water equation
! (SWE).
!
!   \begin{tabular}{|l|l|} \hline \hline
!       c_sw  &   \\ \hline
!       d_sw  & 
!   \end{tabular}
!
! !REVISION HISTORY:
!   01.01.15   Lin        Routines coalesced into this module
!   03.11.19   Sawyer     Merged in CAM changes by Mirin
!   04.10.07   Sawyer     ompinner now from dynamics_vars
!   05.03.25   Todling    shr_kind_r8 can only be referenced once (MIPSpro-7.4.2)
!   05.05.25   Sawyer     Merged CAM and GEOS5 versions (mostly CAM)
!   05.07.26   Worley     Changes for Cray X1
!   05.07.05   Sawyer     Interfaces of c_sw and d_sw simplified with grid
!   05.10.12   Worley     More changes for Cray X1(E), avoiding array segment copying
!   06.01.18   Putman     Allowed Y-dir courant number and mass flux to accumulate
!                         at jlast+1
!   06.09.06   Sawyer     Isolated magic numbers as F90 parameters
!
!EOP

! Magic numbers used in this module

      real(r8), parameter ::  D0_0                    =  0.0_r8
      real(r8), parameter ::  D0_125                  =  0.125_r8
      real(r8), parameter ::  D0_25                   =  0.25_r8
      real(r8), parameter ::  D0_5                    =  0.5_r8
      real(r8), parameter ::  D1_0                    =  1.0_r8
      real(r8), parameter ::  D2_0                    =  2.0_r8
      real(r8), parameter ::  D1E30                   =  1.0e30_r8

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: c_sw --- Solve the SWE on a C grid
!
! !INTERFACE:
 subroutine c_sw(grid,   u,       v,      pt,       delp,               &
                 u2,     v2,                                            &
                 uc,     vc,      ptc,    delpf,    ptk,                &
                 tiny,   iord,    jord)

! Routine for shallow water dynamics on the C-grid

! !USES:

  use tp_core
  use pft_module, only : pft2d

  implicit none

! !INPUT PARAMETERS:
  type (T_FVDYCORE_GRID), intent(in) :: grid
  integer, intent(in):: iord
  integer, intent(in):: jord

  real(r8), intent(in):: u2(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
  real(r8), intent(in):: v2(grid%im,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d)

! Prognostic variables:
  real(r8), intent(in):: u(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_s)
  real(r8), intent(in):: v(grid%im,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d)
  real(r8), intent(in):: pt(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
  real(r8), intent(in):: delp(grid%im,grid%jfirst:grid%jlast)
  real(r8), intent(in):: delpf(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)

  real(r8), intent(in):: tiny

! !INPUT/OUTPUT PARAMETERS:
  real(r8), intent(inout):: uc(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
  real(r8), intent(inout):: vc(grid%im,grid%jfirst-2:grid%jlast+2 )

! !OUTPUT PARAMETERS:
  real(r8), intent(out):: ptc(grid%im,grid%jfirst:grid%jlast)
  real(r8), intent(out):: ptk(grid%im,grid%jfirst:grid%jlast)

! !DESCRIPTION:
!
!   Routine for shallow water dynamics on the C-grid
!
! !REVISION HISTORY:
!   WS   2003.11.19     Merged in CAM changes by Mirin
!   WS   2004.10.07     Added ProTeX documentation
!   WS   2005.07.01     Simplified interface by passing grid
!
!EOP
!-----------------------------------------------------------------------
!BOC


!--------------------------------------------------------------
! Local 
  real(r8) :: zt_c
  real(r8) :: dydt
  real(r8) :: dtdy5
  real(r8) :: rcap

  real(r8), pointer:: sc(:)
  real(r8), pointer:: dc(:,:)

  real(r8), pointer:: cosp(:)
  real(r8), pointer:: acosp(:)
  real(r8), pointer:: cose(:)

  real(r8), pointer:: dxdt(:)
  real(r8), pointer:: dxe(:)
  real(r8), pointer:: rdxe(:)
  real(r8), pointer:: dtdx2(:)
  real(r8), pointer:: dtdx4(:)
  real(r8), pointer:: dtxe5(:)
  real(r8), pointer:: dycp(:)
  real(r8), pointer::  cye(:)

  real(r8), pointer:: fc(:)

  real(r8), pointer:: sinlon(:)
  real(r8), pointer:: coslon(:)
  real(r8), pointer:: sinl5(:)
  real(r8), pointer:: cosl5(:)

    real(r8) :: fx(grid%im,grid%jfirst:grid%jlast)
    real(r8) :: xfx(grid%im,grid%jfirst:grid%jlast)
    real(r8) :: tm2(grid%im,grid%jfirst:grid%jlast)

    real(r8) :: va(grid%im,grid%jfirst-1:grid%jlast)

    real(r8) :: wk4(grid%im+2,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d)
    
    real(r8) :: wk1(grid%im,grid%jfirst-1:grid%jlast+1)

    real(r8) :: cry(grid%im,grid%jfirst-1:grid%jlast+1)
    real(r8) :: fy(grid%im,grid%jfirst-1:grid%jlast+1)
    
    real(r8) :: ymass(grid%im,grid%jfirst: grid%jlast+1) 
    real(r8) :: yfx(grid%im,grid%jfirst: grid%jlast+1)

    real(r8) :: crx(grid%im,grid%jfirst-grid%ng_c:grid%jlast+grid%ng_c)
    real(r8) :: vort_u(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
    real(r8) :: vort(grid%im,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d)

    real(r8) :: fxjv(grid%im,grid%jfirst-1:grid%jn2g0)
    real(r8) :: p1dv(grid%im,grid%jfirst-1:grid%jn2g0)
    real(r8) :: cx1v(grid%im,grid%jfirst-1:grid%jn2g0)

    real(r8) :: qtmp(-grid%im/3:grid%im+grid%im/3)
    real(r8) :: qtmpv(-grid%im/3:grid%im+grid%im/3, grid%jfirst-1:grid%jn2g0)
    real(r8) :: slope(-grid%im/3:grid%im+grid%im/3)
    real(r8) :: al(-grid%im/3:grid%im+grid%im/3)
    real(r8) :: ar(-grid%im/3:grid%im+grid%im/3)
    real(r8) :: a6(-grid%im/3:grid%im+grid%im/3)

    real(r8) :: us, vs, un, vn
    real(r8) :: p1ke, p2ke
    real(r8) :: uanp(grid%im), uasp(grid%im), vanp(grid%im), vasp(grid%im)

    logical :: ffsl(grid%jm)
    logical :: sldv(grid%jfirst-1:grid%jn2g0)

    integer :: i, j, im2
    integer :: js1g1, js2g1, js2gc1, jn2gc, jn1g1, js2g0, js2gc, jn1gc
    integer :: im, jm, jfirst, jlast, jn2g0, ng_s, ng_c, ng_d



!
!   For convenience
!

  im     = grid%im
  jm     = grid%jm
  jfirst = grid%jfirst
  jlast  = grid%jlast

  jn2g0  = grid%jn2g0

  ng_c   = grid%ng_c
  ng_d   = grid%ng_d
  ng_s   = grid%ng_s

  rcap   = grid%rcap

  zt_c   =  grid%zt_c
  dydt   =  grid%dydt
  dtdy5  =  grid%dtdy5

  sc     => grid%sc
  dc     => grid%dc

  cosp   => grid%cosp
  acosp  => grid%acosp
  cose   => grid%cose

  dxdt   => grid%dxdt
  dxe    => grid%dxe
  rdxe   => grid%rdxe
  dtdx2  => grid%dtdx2
  dtdx4  => grid%dtdx4
  dtxe5  => grid%dtxe5
  dycp   => grid%dycp
  cye    => grid%cye
  fc     => grid%fc

  sinlon => grid%sinlon
  coslon => grid%coslon
  sinl5  => grid%sinl5
  cosl5  => grid%cosl5


! Set loop limits

    im2 = im/2

    js2g0  = max(2,jfirst)
    js2gc  = max(2,jfirst-ng_c) ! NG lats on S (starting at 2)
    jn1gc  = min(jm,jlast+ng_c) ! ng_c lats on N (ending at jm)
    js1g1  = max(1,jfirst-1)
    js2g1  = max(2,jfirst-1)
    jn1g1  = min(jm,jlast+1)
    jn2gc  = min(jm-1,jlast+ng_c)   ! NG latitudes on N (ending at jm-1)
    js2gc1 = max(2,jfirst-ng_c+1)   ! NG-1 latitudes on S (starting at 2) 

! KE at poles
    if ( jfirst-ng_d <= 1 ) then
       p1ke = D0_125*(u2(1, 1)**2 + v2(1, 1)**2)
    endif

    if ( jlast+ng_d >= jm ) then
       p2ke = D0_125*(u2(1,jm)**2 + v2(1,jm)**2)
    endif

        if ( jfirst /= 1 ) then
          do i=1,im
            cry(i,jfirst-1) = dtdy5*vc(i,jfirst-1)
          enddo

        endif

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif

        do j=js2g0,jn1g1                     ! ymass needed on NS
          do i=1,im
               cry(i,j) = dtdy5*vc(i,j)
             ymass(i,j) = cry(i,j)*cose(j)
          enddo
        enddo

! New va definition
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
        do j=js2g1,jn2g0                     ! va needed on S (for YCC, iv==1)
          do i=1,im
            va(i,j) = D0_5*(cry(i,j)+cry(i,j+1))
          enddo
        enddo

! SJL: Check if FFSL integer fluxes need to be computed

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
        do j=js2gc,jn2gc                ! ffsl needed on N*sg S*sg
          do i=1,im
            crx(i,j) = uc(i,j)*dtdx2(j)
          enddo
          ffsl(j) = .false.
          if( cosp(j) < zt_c ) then
            do i=1,im
              if( abs(crx(i,j)) > D1_0 ) then
                ffsl(j) = .true. 
#if ( !defined UNICOSMP ) || ( !defined NEC_SX )
                exit
#endif
              endif
            enddo
          endif
        enddo

! 2D transport of polar filtered delp (for computing fluxes!)
! Update is done on the unfiltered delp

   call tp2c( ptk,  va(1,jfirst),  delpf(1,jfirst-ng_c),    &
              crx(1,jfirst-ng_c), cry(1,jfirst),             &
              im, jm, iord, jord, ng_c, xfx,                 &
              yfx, ffsl, rcap, acosp,                        &
              crx(1,jfirst), ymass, cosp,                    &
              0, jfirst, jlast)

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
   do j=js2g0,jn2g0                      ! xfx not ghosted
      if( ffsl(j) ) then
         do i=1,im
           xfx(i,j) = xfx(i,j)/sign(max(abs(crx(i,j)),tiny),crx(i,j))
         enddo
      endif
   enddo

! pt-advection using pre-computed mass fluxes
! use tm2 below as the storage for pt increment
! WS 99.09.20 : pt, crx need on N*ng S*ng, yfx on N

    call tp2c(tm2 ,va(1,jfirst), pt(1,jfirst-ng_c),       &
              crx(1,jfirst-ng_c), cry(1,jfirst),          &
              im, jm,  iord, jord, ng_c, fx,              &
              fy(1,jfirst), ffsl, rcap, acosp,            &
              xfx, yfx, cosp, 1, jfirst, jlast)

! use wk4, crx as work arrays
     call pft2d(ptk(1,js2g0), sc,   &
                dc, im, jn2g0-js2g0+1,  &
                wk4, crx )
     call pft2d(tm2(1,js2g0), sc,   &
                dc, im, jn2g0-js2g0+1,  &
                wk4, crx )

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
    do j=jfirst,jlast
       do i=1,im
          ptk(i,j) = delp(i,j) + ptk(i,j)
          ptc(i,j) = (pt(i,j)*delp(i,j) + tm2(i,j))/ptk(i,j)
       enddo
    enddo

!------------------
! Momentum equation
!------------------

     call ycc(im, jm, fy, vc(1,jfirst-2), va(1,jfirst-1),   &
              va(1,jfirst-1), jord, 1, jfirst, jlast)

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
     do j=js2g1,jn2g0

          do i=1,im
            cx1v(i,j) = dtdx4(j)*u2(i,j)
          enddo

          sldv(j) = .false.
          if( cosp(j) < zt_c ) then
            do i=1,im
              if( abs(cx1v(i,j)) > D1_0 ) then
                sldv(j) = .true. 
#if ( !defined UNICOSMP ) || ( !defined NEC_SX )
                exit
#endif
              endif
            enddo
          endif

          p1dv(im,j) = uc(1,j)
          do i=1,im-1
            p1dv(i,j) = uc(i+1,j)
          enddo

     enddo

     call xtpv(im, sldv, fxjv, p1dv, cx1v, iord, cx1v,        &
              cosp, 0, slope, qtmpv, al, ar, a6,              &
              jfirst, jlast, js2g1, jn2g0, jm,                &
              jfirst-1, jn2g0, jfirst-1, jn2g0,               &
              jfirst-1, jn2g0, jfirst-1, jn2g0,               &
              jfirst-1, jn2g0, jfirst-1, jn2g0)

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
     do j=js2g1,jn2g0
        do i=1,im
          wk1(i,j) = dxdt(j)*fxjv(i,j) + dydt*fy(i,j)
       enddo
     enddo

     if ( jfirst == 1 ) then
          do i=1,im
            wk1(i,1) = p1ke
          enddo
     endif

     if ( jlast == jm ) then
          do i=1,im
            wk1(i,jm) = p2ke
          enddo
     endif

! crx redefined
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
     do j=js2gc1,jn1gc
            crx(1,j) = dtxe5(j)*u(im,j)
          do i=2,im
            crx(i,j) = dtxe5(j)*u(i-1,j)
          enddo
     enddo

     if ( jfirst /=1 ) then 
          do i=1,im
             cry(i,jfirst-1) = dtdy5*v(i,jfirst-1)
          enddo
     endif

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
     do j=jfirst,jlast
        do i=1,im
             cry(i,j) = dtdy5*v(i,j)
           ymass(i,j) = cry(i,j)*cosp(j)       ! ymass actually unghosted
        enddo
     enddo

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
     do j=js2g0,jlast
          do i=1,im
            tm2(i,j) = D0_5*(cry(i,j)+cry(i,j-1)) ! cry ghosted on S 
          enddo
     enddo

!    Compute absolute vorticity on the C-grid.

     if ( jfirst-ng_d <= 1 ) then
          do i=1,im
            vort_u(i,1) = D0_0
          enddo
     endif

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
     do j=js2gc,jn2gc
         do i=1,im
            vort_u(i,j) = uc(i,j)*cosp(j)
         enddo
     enddo

     if ( jlast+ng_d >= jm ) then
          do i=1,im
            vort_u(i,jm) = D0_0
          enddo
     endif

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
     do j=js2gc1,jn1gc
! The computed absolute vorticity on C-Grid is assigned to vort
          vort(1,j) = fc(j) + (vort_u(1,j-1)-vort_u(1,j))*cye(j) +     &
                    (vc(1,j) - vc(im,j))*rdxe(j)

          do i=2,im
             vort(i,j) = fc(j) + (vort_u(i,j-1)-vort_u(i,j))*cye(j) +  &
                       (vc(i,j) - vc(i-1,j))*rdxe(j)
          enddo
     enddo

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
     do j=js2gc1,jn1gc          ! ffsl needed on N*ng S*(ng-1)
          ffsl(j) = .false.
          if( cose(j) < zt_c ) then
            do i=1,im
              if( abs(crx(i,j)) > D1_0 ) then
                ffsl(j) = .true. 
#if ( !defined UNICOSMP ) || ( !defined NEC_SX )
                exit
#endif
              endif
            enddo
          endif
     enddo

   call tpcc( tm2, ymass, vort(1,jfirst-ng_d), crx(1,jfirst-ng_c),  &
              cry(1,jfirst), im, jm, ng_c, ng_d,                  &
              iord, jord, fx, fy(1,jfirst), ffsl, cose,           &
              jfirst, jlast, slope, qtmp, al, ar, a6 )

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
   do j=js2g0,jn2g0
         uc(1,j) = uc(1,j) + dtdx2(j)*(wk1(im,j)-wk1(1,j)) + dycp(j)*fy(1,j)
      do i=2,im
         uc(i,j) = uc(i,j) + dtdx2(j)*(wk1(i-1,j)-wk1(i,j)) + dycp(j)*fy(i,j)
      enddo
   enddo
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
   do j=js2g0,jlast
        do i=1,im-1
           vc(i,j) = vc(i,j) + dtdy5*(wk1(i,j-1)-wk1(i,j))-dxe(j)*fx(i+1,j)
        enddo
           vc(im,j) = vc(im,j) + dtdy5*(wk1(im,j-1)-wk1(im,j))-dxe(j)*fx(1,j)
   enddo
!EOC
 end subroutine c_sw
!--------------------------------------------------------------------------



!-----------------------------------------------------------------------
!BOP
! !IROUTINE: d_sw --- Solve the SWE on a D grid
!
! !INTERFACE:
 subroutine d_sw( grid,  u,      v,        uc,     vc,               &   
                  pt,    delp,   delpf,    cx3,    cy3,              &
                  mfx,   mfy,    cdx,      cdy,          &
                  cdxde, cdxdp,  cdyde,   cdydp,                     & !ldel2 variables
                  cdxdiv,  cdydiv, cdx4,  cdy4,  cdtau4,  &
                  ldiv2, ldiv4, ldel2, &
                  iord,  jord,  tiny )  
!--------------------------------------------------------------------------
! Routine for shallow water dynamics on the D-grid

! !USES:

  use tp_core
  use pft_module, only : pft2d

  implicit none

! !INPUT PARAMETERS:
  type (T_FVDYCORE_GRID), intent(in) :: grid
  integer, intent(in):: iord
  integer, intent(in):: jord
  logical, intent(in) :: ldiv2,ldiv4,ldel2   !damping options

! Prognostic variables:
  real(r8), intent(inout):: u(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_s)
  real(r8), intent(inout):: v(grid%im,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d)
! Delta pressure
  real(r8), intent(inout):: delp(grid%im,grid%jfirst:grid%jlast)
! Potential temperature
  real(r8), intent(inout):: pt(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)

  real(r8), intent(inout):: delpf(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)

  real(r8), intent(in):: cdx   (grid%js2g0:grid%jn1g1)
  real(r8), intent(in):: cdy   (grid%js2g0:grid%jn1g1)
  !
  ! variables for div4 and del2 damping
  !
  real(r8), intent(in):: cdx4  (grid%js2g0:grid%jn1g1) 
  real(r8), intent(in):: cdy4  (grid%js2g0:grid%jn1g1) 
  real(r8), intent(in):: cdtau4(grid%js2g0:grid%jn1g1) 
  real(r8), intent(in):: cdxde (grid%js2g0:grid%jn1g1) 
  real(r8), intent(in):: cdxdp (grid%js2g0:grid%jn1g1) 
  real(r8), intent(in):: cdydp (grid%js2g0:grid%jn1g1) 
  real(r8), intent(in):: cdyde (grid%js2g0:grid%jn1g1) 
  real(r8), intent(in):: cdxdiv(grid%jm)          
  real(r8), intent(in):: cdydiv(grid%jm)          

  real(r8), intent(in):: tiny

! !INPUT/OUTPUT PARAMETERS:
  real(r8), intent(inout):: uc(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
  real(r8), intent(inout):: vc(grid%im,grid%jfirst-2   :grid%jlast+2 )
  real(r8), intent(inout):: cx3(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)! Accumulated Courant number in X
  real(r8), intent(inout):: cy3(grid%im,grid%jfirst:grid%jlast+1)        ! Accumulated Courant number in Y
  real(r8), intent(inout):: mfx(grid%im,grid%jfirst:grid%jlast)          ! Mass flux in X  (unghosted)
  real(r8), intent(inout):: mfy(grid%im,grid%jfirst:grid%jlast+1)        ! Mass flux in Y

! !DESCRIPTION:
!
!   Routine for shallow water dynamics on the D-grid
!
! !REVISION HISTORY:
!   WS   2003.11.19     Merged in CAM changes by Mirin
!   WS   2004.10.07     Added ProTeX documentation
!   WS   2005.07.05     Simplified interface using grid
!
!EOP
!-----------------------------------------------------------------------
!BOC


! Local
  integer :: im
  integer :: jm
  integer :: jfirst
  integer :: jlast
  integer :: js2g0
  integer :: jn1g1
  integer :: ng_d
  integer :: ng_s
  integer :: nq

  real(r8) :: zt_d
  real(r8) :: tdy5
  real(r8) :: rdy
  real(r8) :: dtdy
  real(r8) :: dtdy5
  real(r8) :: rcap

  real(r8), pointer:: sc(:)
  real(r8), pointer:: dc(:,:)
  real(r8), pointer:: se(:)
  real(r8), pointer:: de(:,:)

  real(r8), pointer :: cosp(:)
  real(r8), pointer :: acosp(:)
  real(r8), pointer :: cose(:)

  real(r8), pointer :: sinlon(:)
  real(r8), pointer :: coslon(:)
  real(r8), pointer :: sinl5(:)
  real(r8), pointer :: cosl5(:)

  real(r8), pointer :: dtdx(:)
  real(r8), pointer :: dtdxe(:)
  real(r8), pointer :: dx(:)
  real(r8), pointer :: rdx(:)
  real(r8), pointer :: cy(:)
  real(r8), pointer :: dyce(:)
  real(r8), pointer :: dtxe5(:)
  real(r8), pointer :: txe5(:)

  real(r8), pointer :: f0(:)
  
  real(r8)   fx(grid%im,grid%jfirst:grid%jlast)
  real(r8)  xfx(grid%im,grid%jfirst:grid%jlast)
  
  !for del2 damping
  real(r8) :: wku(grid%im,grid%jfirst-1:grid%jlast+1) 
  real(r8) :: wkv(grid%im,grid%jfirst-1:grid%jlast+1) 
 
  !for div4 damping
  real(r8) ::   wkdiv4(grid%im+2,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_s) 
  real(r8) ::  wk2div4(grid%im+1,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_s) 
  
  real(r8)  wk1(grid%im,grid%jfirst-1:grid%jlast+1)
  real(r8)  cry(grid%im,grid%jfirst-1:grid%jlast+1)
  real(r8)   fy(grid%im,grid%jfirst-2:grid%jlast+2)!halo changed for div4
  
  real(r8) ymass(grid%im,grid%jfirst: grid%jlast+1) 
  real(r8)   yfx(grid%im,grid%jfirst: grid%jlast+1)
  
  real(r8)   va(grid%im,grid%jfirst-1:grid%jlast)
  real(r8)   ub(grid%im,grid%jfirst:  grid%jlast+1)
  
  real(r8)  crx(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
#if defined(FILTER_MASS_FLUXES)
  real(r8)   u2(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
  real(r8)   v2(grid%im+2,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d)
#endif
  
  real(r8) fxjv(grid%im,grid%jfirst-1:grid%jn1g1)
  real(r8) qtmpv(-grid%im/3:grid%im+grid%im/3, grid%jfirst-1:grid%jn1g1)
  real(r8) slope(-grid%im/3:grid%im+grid%im/3)
  real(r8) al(-grid%im/3:grid%im+grid%im/3)
  real(r8) ar(-grid%im/3:grid%im+grid%im/3)
  real(r8) a6(-grid%im/3:grid%im+grid%im/3)
  
  real(r8)  c1, c2
  real(r8) uanp(grid%im), uasp(grid%im), vanp(grid%im), vasp(grid%im)
  real(r8) un, vn, us, vs, r2im
  
  real(r8)  div (grid%im,grid%jfirst-1:grid%jlast+2) !for div4 damping
  real(r8)  div4(grid%im,grid%jfirst-1:grid%jlast+1) !for div4 damping
  
  logical ffsl(grid%jm)
  logical sldv(grid%jfirst-1:grid%jn1g1)
  
  real(r8):: deldiv     !for div4
  
  
  integer i, j
  integer js2gd, jn2g0, jn2g1, jn2gd, jn1gd
  integer jn2g2 !for extra halo for div4 
  integer js2gs, jn2gs, jn1gs
  integer im2
  
!
! For convenience
!
  nq     =  grid%nq

  im     =  grid%im
  jm     =  grid%jm
  jfirst =  grid%jfirst
  jlast  =  grid%jlast
  ng_d   =  grid%ng_d
  ng_s   =  grid%ng_s
  js2g0  =  grid%js2g0
  jn1g1  =  grid%jn1g1

  rcap   =  grid%rcap
  zt_d   =  grid%zt_d
  tdy5   =  grid%tdy5
  rdy    =  grid%rdy
  dtdy   =  grid%dtdy
  dtdy5  =  grid%dtdy5

  sc     => grid%sc
  dc     => grid%dc
  se     => grid%se
  de     => grid%de

  cosp   => grid%cosp
  acosp  => grid%acosp
  cose   => grid%cose

  sinlon => grid%sinlon
  coslon => grid%coslon
  sinl5  => grid%sinl5
  cosl5  => grid%cosl5

  dtdx   => grid%dtdx
  dtdxe  => grid%dtdxe
  dx     => grid%dx
  rdx    => grid%rdx
  cy     => grid%cy
  dyce   => grid%dyce
  dtxe5  => grid%dtxe5
  txe5   => grid%txe5

  f0     => grid%f0

! Set loop limits

  jn2g0 = min(jm-1,jlast)
  jn2g1 = min(jm-1,jlast+1)
  jn2g2 = min(jm-1,jlast+2) 

  js2gd = max(2,jfirst-ng_d)     ! NG latitudes on S (starting at 1)
  jn2gd = min(jm-1,jlast+ng_d)   ! NG latitudes on S (ending at jm-1)
  jn1gd = min(jm,jlast+ng_d)     ! NG latitudes on N (ending at jm)
  js2gs = max(2,jfirst-ng_s)
  jn2gs = min(jm-1,jlast+ng_s)
  jn1gs = min(jm,jlast+ng_s)     ! NG latitudes on N (ending at jm)

! Get C-grid U-wind at poles.
  im2 = im/2
  r2im = 0.5_r16/real(im,r16)

  if ( jfirst <= 1 ) then
!
! Treat SP
!
     do i=1,im-1
        uasp(i) = uc(i,2) + uc(i+1,2)
        vasp(i) = vc(i,2) + vc(i,3)
     enddo
     uasp(im) = uc(im,2) + uc(1,2)
     vasp(im) = vc(im,2) + vc(im,3)

! Projection at SP

     us = D0_0
     vs = D0_0
     do i=1,im2
        us = us + (uasp(i+im2)-uasp(i))*sinlon(i)     &
                + (vasp(i)-vasp(i+im2))*coslon(i)
        vs = vs + (uasp(i+im2)-uasp(i))*coslon(i)     &
                + (vasp(i+im2)-vasp(i))*sinlon(i)
     enddo
     us = us*r2im
     vs = vs*r2im

! get U-wind at SP

     do i=1,im2
        uc(i,  1) = -us*sinl5(i) - vs*cosl5(i)
        uc(i+im2,  1) = -uc(i,  1)
     enddo

  endif

  if ( jlast >= jm ) then
!
! Treat NP
!
     do i=1,im-1
        uanp(i) = uc(i,jm-1) + uc(i+1,jm-1)
        vanp(i) = vc(i,jm-1) + vc(i,jm)
     enddo
     uanp(im) = uc(im,jm-1) + uc(1,jm-1)
     vanp(im) = vc(im,jm-1) + vc(im,jm)

! Projection at NP

     un = D0_0
     vn = D0_0
     do i=1,im2
        un = un + (uanp(i+im2)-uanp(i))*sinlon(i)  &
                + (vanp(i+im2)-vanp(i))*coslon(i)
        vn = vn + (uanp(i)-uanp(i+im2))*coslon(i)  &
                + (vanp(i+im2)-vanp(i))*sinlon(i)
     enddo
     un = un*r2im
     vn = vn*r2im

! get U-wind at NP

     do i=1,im2
        uc(i,jm) = -un*sinl5(i) + vn*cosl5(i)
        uc(i+im2,jm) = -uc(i,jm)
     enddo

  endif

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
  do j=js2gd,jn2gd                     ! crx needed on N*ng S*ng
     do i=1,im
        crx(i,j) = dtdx(j)*uc(i,j)
     enddo
  enddo

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
  do j=js2gd,jn2gd                ! ffsl needed on N*ng S*ng
     ffsl(j) = .false.
     if( cosp(j) < zt_d ) then
         do i=1,im
            if( abs(crx(i,j)) > D1_0 ) then
               ffsl(j) = .true. 
#if ( !defined UNICOSMP ) || ( !defined NEC_SX )
               exit
#endif
            endif
         enddo
      endif
  enddo

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
  do j=js2g0,jn1g1                       ! cry, ymass needed on N
     do i=1,im
        cry(i,j) = dtdy*vc(i,j)
        ymass(i,j) = cry(i,j)*cose(j)
     enddo
  enddo

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
  do j=js2g0,jn2g0                         ! No ghosting
     do i=1,im
        if( cry(i,j)*cry(i,j+1) > D0_0 ) then
           if( cry(i,j) > D0_0 ) then
              va(i,j) = cry(i,j)
           else
              va(i,j) = cry(i,j+1)         ! cry ghosted on N
           endif
        else
           va(i,j) = D0_0
        endif
     enddo
  enddo

! transport polar filtered delp
      call tp2c(ub(1,jfirst), va(1,jfirst), delpf(1,jfirst-ng_d),   &
                crx(1,jfirst-ng_d),cry(1,jfirst),im,jm,iord,jord,   &
                ng_d, xfx, yfx, ffsl,                               &
                rcap, acosp,crx(1,jfirst), ymass,                   &
                cosp, 0, jfirst, jlast)

#if defined(FILTER_MASS_FLUXES)
   call pft2d( xfx(1,js2g0), sc, dc, im, jn2g0-js2g0+1, &
                    v2, u2 )
   call pft2d( yfx(1,js2g0), se, de, im, jn1g1-js2g0+1, &
                    v2, u2 )
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
   do j=js2g0,jn2g0
      do i=1,im-1
         ub(i,j) = xfx(i,j) - xfx(i+1,j) + (yfx(i,j)-yfx(i,j+1))*acosp(j)
      enddo
      ub(im,j) = xfx(im,j) - xfx(1,j) + (yfx(im,j)-yfx(im,j+1))*acosp(j)
   enddo
#endif

! <<< Save necessary data for large time step tracer transport >>>
      if( nq > 0 ) then
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
          do j=js2g0,jn2g0                       ! No ghosting needed
            do i=1,im
              cx3(i,j) = cx3(i,j) + crx(i,j)
              mfx(i,j) = mfx(i,j) + xfx(i,j)
            enddo
          enddo

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
          do j=js2g0,jlast                      ! No ghosting needed
            do i=1,im
              cy3(i,j) = cy3(i,j) + cry(i,j)
              mfy(i,j) = mfy(i,j) + yfx(i,j)
            enddo
          enddo
      endif
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
     do j=js2g0,jn2g0                         ! No ghosting needed
        if( ffsl(j) ) then
          do i=1,im
             xfx(i,j) = xfx(i,j)/sign(max(abs(crx(i,j)),tiny),crx(i,j))
          enddo
        endif
     enddo

! Update delp
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
        do j=jfirst,jlast
          do i=1,im
! SAVE old delp: pressure thickness ~ "air density"
            wk1(i,j) = delp(i,j)
            delp(i,j) = wk1(i,j) + ub(i,j)
          enddo
        enddo

! pt Advection
  call tp2c(ub(1,jfirst),va(1,jfirst),pt(1,jfirst-ng_d),    &
            crx(1,jfirst-ng_d),cry(1,jfirst),               &
            im,jm,iord,jord,ng_d,fx,fy(1,jfirst),           &
            ffsl, rcap, acosp,                              &
            xfx, yfx(1,jfirst), cosp, 1, jfirst,jlast)

! Update pt.
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=jfirst,jlast
         do i=1,im
            pt(i,j) = (pt(i,j)*wk1(i,j)+ub(i,j)) / delp(i,j)
         enddo
      enddo

! Compute upwind biased kinetic energy at the four cell corners

! Start using ub as v (CFL) on B-grid (cell corners)
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js2g0,jn1g1                          ! ub needed on N
           ub(1,j) = dtdy5*(vc(1,j) + vc(im,j))  
         do i=2,im
            ub(i,j) = dtdy5*(vc(i,j) + vc(i-1,j))
         enddo
      enddo

      call ytp(im, jm, fy(1,jfirst), v(1,jfirst-ng_d), ub(1,jfirst),  &
               ub(1,jfirst), ng_d, jord, 1, jfirst, jlast)
! End using ub as v (CFL) on B-grid

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
   do j=js2g0,jn1g1                 ! ub needed on N
       do i=1,im                
          ub(i,j) = dtxe5(j)*(uc(i,j) + uc(i,j-1))
!                        uc will be used as wrok array after this point
       enddo
   enddo

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
  do j=js2g0,jn1g1                       ! wk1 needed on N
          sldv(j) = .false.
          if( cose(j) < zt_d ) then
            do i=1,im
              if( abs(ub(i,j)) > D1_0 ) then    ! ub ghosted on N
                sldv(j) = .true. 
#if ( !defined UNICOSMP ) || ( !defined NEC_SX )
                exit
#endif
              endif
            enddo
          endif
  enddo

  call xtpv(im,  sldv, fxjv, u, ub, iord, ub, cose,       &
           0, slope, qtmpv, al, ar, a6,                   &
           jfirst, jlast, js2g0, jn1g1, jm,               &
           jfirst-1, jn1g1, jfirst-1, jn1g1,              &
           jfirst-ng_d, jlast+ng_s, jfirst, jlast+1,      &
           jfirst, jlast+1, jfirst-1, jn1g1) 
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
  do j=js2g0,jn1g1                       ! wk1 needed on N
     do i=1,im
        wk1(i,j) =  txe5(j)*fxjv(i,j) + tdy5*fy(i,j)  ! fy ghosted on N
     enddo
  enddo

! Add divergence damping to vector invariant form of the momentum eqn
! (absolute vorticity is damped by ffsl scheme, therefore divergence damping
! provides more consistent dissipation to divergent part of the flow)

!--------------------------
! Perform divergence damping 
!--------------------------

  if (ldiv2) then
     !
     ! standard div2 damping
     !
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
     do j=max(2,jfirst-1), jn2g1                   ! fy need on NS (below)
        do i=1,im
           !
           ! cosp is cosine(theta) at cell center disctretized from the identity 
           !
           !   cos(theta) = d(sin(theta))/d(theta)
           !
           ! as
           !
           !   cosp = (sine(j+1)-sine(j))/dp where dp = pi/(jm-1)
           !
           fy(i,j) = v(i,j)*cosp(j)      ! v ghosted on NS at least
        enddo
     enddo

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif   
     do j=js2g0,jn1g1
        ! i=1
        uc(1,j) = u(im,j) - u(1,j)    ! u ghosted on N at least
        do i=2,im
           uc(i,j) = u(i-1,j) - u(i,j)
        enddo
     enddo
     
     if ( jfirst == 1 ) then
        ! j=2
        do i=1,im
           wk1(i,2) = wk1(i,2) - cdy(2)*fy(i, 2) + cdx(2)*uc(i,2)
        enddo
     endif
     
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif     
     do j=max(3,jfirst),jn2g1         ! wk1 needed on N (after TP2D)
        do i=1,im
           wk1(i,j) = wk1(i,j) + cdy(j)*(fy(i,j-1) - fy(i,j))  &
                + cdx(j)*uc(i,j)
        enddo
     enddo
     
     if ( jlast == jm ) then
        do i=1,im
           wk1(i,jm) = wk1(i,jm) + cdy(jm)*fy(i,jm-1) + cdx(jm)*uc(i,jm)
        enddo
     endif
  end if
  !
   !
   ! js2gd = max(2,jfirst-ng_d)     ! NG latitudes on S (starting at 1)
   ! jn2gd = min(jm-1,jlast+ng_d)   ! NG latitudes on S (ending at jm-1)
   ! jn1gd = min(jm,jlast+ng_d)     ! NG latitudes on N (ending at jm)
   ! js2gs = max(2,jfirst-ng_s)
   ! jn2gs = min(jm-1,jlast+ng_s)
   ! jn1gs = min(jm,jlast+ng_s)     ! NG latitudes on N (ending at jm)
  if (ldiv4) then
     !
     ! filter velocity components for stability
     !
     call pft2d(u(1,js2gd), grid%sediv4, grid%dediv4, im, jn1gs-js2gd+1, &
          wkdiv4, wk2div4 )
    
     call pft2d(v(1,js2gs), grid%scdiv4, grid%dcdiv4, im, jn2gd-js2gs+1, &
          wkdiv4, wk2div4 )


    !**************************************************************************
    !
    ! div4 damping
    !
    !**************************************************************************

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
    do j=max(2,jfirst-2), min(jm-1,grid%jlast+2)                   ! fy need on NS (below)
      do i=1,im
        fy(i,j) = v(i,j)*cosp(j)      ! v ghosted on NS at least
      enddo
    enddo
    
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif   
    do j=max(2,jfirst-1),min(jm,grid%jlast+2)
      ! i=1
      uc(1,j) = u(im,j) - u(1,j)    ! u ghosted on N at least
      do i=2,im
        uc(i,j) = u(i-1,j) - u(i,j)
      enddo
    enddo
    !
    ! compute divergence
    !
    if ( jfirst == 1 ) then
      ! j=2
      do i=1,im
        div(i,2) = - cdydiv(2)*fy(i, 2) + cdxdiv(2)*uc(i,2)
      enddo
    endif
    
#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif     
   do j=max(3,jfirst-1),min(jm-1,grid%jlast+2)!jn2g2         ! wk1 needed on N (after TP2D)
      do i=1,im
        div(i,j) = cdydiv(j)*(fy(i,j-1) - fy(i,j)) + cdxdiv(j)*uc(i,j)
      enddo
    enddo
    
    if ( jlast == jm ) then
      do i=1,im
        div(i,jm) = cdydiv(jm)*fy(i,jm-1) + cdxdiv(jm)*uc(i,jm)
      enddo
    endif
    
    if ( jfirst == 1 ) then
          ! j=2
      i=1
      j=2
      deldiv = cdx4(j) * (div(i+1,j  )-D2_0*div(i,j)+div(im ,j  ))+&
               cdy4(j) * (cosp(j)*(div(i,j+1)-div(i,j)))
      wk1(i,j) = wk1(i,j) +cdtau4(j)*deldiv
      do i=2,im-1
        deldiv = cdx4(j) * (div(i+1,j  )-D2_0*div(i,j)+div(i-1,j  ))+&
                 cdy4(j) * (cosp(j  )*(div(i  ,j+1)-div(i         ,j)))
        wk1(i,j) = wk1(i,j) + cdtau4(j)*deldiv
      enddo
      i=im
      deldiv = cdx4(j) * (div(1 ,j  )-D2_0*div(i,j)+div(i-1,j  ))+&
               cdy4(j) * (cosp(j  )*(div(i,j+1)-div(i,j)))
      wk1(i,j) = wk1(i,j) + cdtau4(j)*deldiv
    endif
    
    do j=max(3,jfirst),jn2g1         ! wk1 needed on N (after TP2D)
      i=1
      deldiv = cdx4(j) * (div(i+1,j  )-D2_0*div(i,j)+div(im ,j  ))+&
               cdy4(j) * ( &
                           cosp(j  )*(div(i  ,j+1)-div(i,j  ))-&
                           cosp(j-1)*(div(i  ,j  )-div(i,j-1)))
      wk1(i,j) = wk1(i,j) +cdtau4(j)*deldiv
      do i=2,im-1
        deldiv = cdx4(j) * (div(i+1,j  )-D2_0*div(i,j)+div(i-1,j  ))+&
                 cdy4(j) * (  &
                             cosp(j  )*(div(i  ,j+1)-div(i         ,j  ))-&
                             cosp(j-1)*(div(i  ,j  )-div(i         ,j-1)))
        wk1(i,j)   = wk1(i,j) + cdtau4(j)*deldiv
      enddo
      i=im
      deldiv = cdx4(j) * (div(1 ,j  )-D2_0*div(i,j)+div(i-1,j  ))+&
               cdy4(j) * ( &
                             cosp(j  )*(div(i  ,j+1)-div(i,j  ))-&
                             cosp(j-1)*(div(i  ,j  )-div(i,j-1)))
      wk1(i,j) = wk1(i,j) + cdtau4(j)*deldiv
    enddo

    if ( jlast == jm ) then
      i=1
      j = jm
      deldiv = cdx4(j) * (div(i+1,j  )-D2_0*div(i,j)+div(im,j  ))+&
               cdy4(j) * (-cosp(j-1)*(div(i,j)-div(i,j-1)))
                       
      wk1(i,j) = wk1(i,j) +cdtau4(j)*deldiv

      do i=2,im-1
        deldiv = cdx4(j) * (div(i+1,j  )-D2_0*div(i,j)+div(i-1,j  ))+&
                 cdy4(j) * (-cosp(j-1)*(div(i,j)-div(i,j-1)))
        wk1(i,j) = wk1(i,j) + cdtau4(j)*deldiv
      enddo
      i=im
      j=jm
      deldiv = cdx4(j) * (div(1,j  )-D2_0*div(i,j)+div(i-1,j  ))+&
               cdy4(j) * (-cosp(j-1)*(div(i,j)-div(i,j-1)))
      wk1(i,j) = wk1(i,j) +cdtau4(j)*deldiv
    endif
  end if

  wku(:,:) = D0_0
  wkv(:,:) = D0_0
  if (ldel2) then
    !**************************************************************************
    !
    ! regular del2 (Laplacian) damping
    !
    !**************************************************************************    
    if ( jfirst == 1 ) then
       ! j=2
       i=1
       j=2

       wku(i,j) = cdxde(j)* (u(i+1,j  )-D2_0*u(i,j)+u(im ,j  ))+&
                  cdyde(j)* (cosp(j  )*(u(i  ,j+1)-u(i         ,j)))
       wkv(i,j) = cdxdp(j) * (v(i+1,j  )-D2_0*v(i,j)+v(im,j  ))+&
                  cdydp(j) * ( &
                           cose(j+1)*(v(i  ,j+1)-v(i,j  ))-&
                           cose(j  )*(v(i  ,j  )        ))
       !line above: there is no v(i,j-1) since it is on the pole
       do i=2,im-1
          wku(i,j) = cdxde(j) * (u(i+1,j  )-D2_0*u(i,j)+u(i-1,j  ))+&
                     cdyde(j) * (cosp(j  )*(u(i  ,j+1)-u(i         ,j)))
          wkv(i,j) = cdxdp(j) * (v(i+1,j  )-D2_0*v(i,j)+v(i-1,j  ))+&
                     cdydp(j) * ( &
                                 cose(j+1)*(v(i  ,j+1)-v(i,j  ))-&
                                 cose(j  )*(v(i  ,j  )        ))        
      enddo
      i=im
      wku(i,j) = cdxde(j) * (u(1 ,j  )-D2_0*u(i,j)+u(i-1,j  ))+&
                 cdyde(j) * (cosp(j  )*(u(i  ,j+1)-u(i         ,j)))                             
      wkv(i,j) = cdxdp(j) * (v(1,j  )-D2_0*v(i,j)+v(i-1 ,j  ))+&
                 cdydp(j) * ( &
                           cose(j+1)*(v(i  ,j+1)-v(i,j  ))-&
                           cose(j  )*(v(i  ,j  )        ))
    endif
    
    do j=max(3,jfirst),jn2g1         ! wk1 needed on N (after TP2D)
      i=1
      wku(i,j) = cdxde(j) * (u(i+1,j  )-D2_0*u(i,j)+u(im ,j  ))+&
                 cdyde(j) * ( &
                           cosp(j  )*(u(i  ,j+1)-u(i,j  ))-&
                           cosp(j-1)*(u(i  ,j  )-u(i,j-1)))

      wkv(i,j) = cdxdp(j) * (v(i+1,j  )-D2_0*v(i,j)+v(im ,j  ))+&
                 cdydp(j) * ( &
                           cose(j+1)*(v(i  ,j+1)-v(i,j  ))-&
                           cose(j  )*(v(i  ,j  )-v(i,j-1)))
      do i=2,im-1
        wku(i,j) = cdxde(j) * (u(i+1,j  )-D2_0*u(i,j)+u(i-1,j  ))+&
                   cdyde(j) * (  &
                             cosp(j  )*(u(i  ,j+1)-u(i         ,j  ))-&
                             cosp(j-1)*(u(i  ,j  )-u(i         ,j-1)))

        wkv(i,j) = cdxdp(j) * (v(i+1,j  )-D2_0*v(i,j)+v(i-1,j  ))+&
                   cdydp(j) * (  &
                             cose(j+1)*(v(i  ,j+1)-v(i         ,j  ))-&
                             cose(j  )*(v(i  ,j  )-v(i         ,j-1)))
      enddo
      i=im
      wku(i,j) = cdxde(j) * (u(1 ,j  )-D2_0*u(i,j)+u(i-1,j  ))+&
                 cdyde(j) * ( &
                             cosp(j  )*(u(i  ,j+1)-u(i,j  ))-&
                             cosp(j-1)*(u(i  ,j  )-u(i,j-1)))

      wkv(i,j) = cdxdp(j) * (v(1 ,j  )-D2_0*v(i,j)+v(i-1,j  ))+&
                 cdydp(j) * ( &
                             cose(j+1)*(v(i  ,j+1)-v(i,j  ))-&
                             cose(j  )*(v(i  ,j  )-v(i,j-1)))
    enddo

    if ( jlast == jm ) then
      i=1
      j = jm
      wku(i,jm) = cdxde(j) * (u(i+1,j  )-D2_0*u(i,j)+u(im,j  ))+&
                  cdyde(j) * (-cosp(j-1)*(u(i,j)-u(i,j-1)))
      do i=2,im-1
         wku(i,jm) = cdxde(j) * (u(i+1,j)-D2_0*u(i,j)+u(i-1,j))+&
                     cdyde(j) * (-cosp(j-1)*(u(i,j)-u(i,j-1)))
      enddo
      i=im
      j=jm
      wku(i,jm) = cdxde(j) * (u(1,j)-D2_0*u(i,j)+u(i-1,j))+&
                  cdyde(j) * (-cosp(j-1)*(u(i,j)-u(i,j-1)))
   endif
end if

!------------------------------------
! End divergence damping computation
!------------------------------------


! Compute Vorticity on the D grid
! delpf used as work array

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js2gd,jn1gd
         do i=1,im
            delpf(i,j) = u(i,j)*cose(j)   ! u ghosted on N*ng S*ng
         enddo
      enddo


      if ( jfirst-ng_d <= 1 ) then
          c1 = D0_0
          do i=1,im
            c1 = c1 + delpf(i,2)
          end do
          c1 = -c1*rdy*rcap

          do i=1,im
            uc(i,1) = c1
          enddo
      endif

      if ( jlast+ng_d >= jm ) then
          c2 = D0_0
          do i=1,im
            c2 = c2 + delpf(i,jm)
          end do
          c2 = c2*rdy*rcap

          do i=1,im
            uc(i,jm) = c2
          enddo
      else

! This is an attempt to avoid ghosting u on N*(ng+1)
          do i=1,im
! DEBUG
!            uc(i,jn2gd) = 0.0
! testing
             uc(i,jn2gd) = D1E30
          enddo
      endif

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js2gd, min(jm-1,jlast+ng_d-1)
          do i=1,im-1
             uc(i,j) = ( delpf(i,j) - delpf(i,j+1)) * cy(j)  +         &
                        (v(i+1,j) - v(i,j))    * rdx(j)
          enddo
            uc(im,j) = (delpf(im,j) - delpf(im,j+1)) *  cy(j) +        &
                       (v(1,j) - v(im,j)) * rdx(j)
      enddo

! uc is relative vorticity at this point

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=max(1,jfirst-ng_d), jn1gd
          do i=1,im
             uc(i,j) = uc(i,j) + f0(j)
! uc is absolute vorticity
          enddo
      enddo

      call tp2d(va(1,jfirst), uc(1,jfirst-ng_d), crx(1,jfirst-ng_d),  &
                cry(1,jfirst), im, jm, iord, jord, ng_d, fx,           &
                fy(1,jfirst), ffsl, crx(1,jfirst),                     &
                ymass, cosp, 0, jfirst, jlast)

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js2g0,jlast
          do i=1,im-1
            uc(i,j) = dtdxe(j)*(wk1(i,j)-wk1(i+1,j)) + dyce(j)*fy(i,j)+wku(i,j)
          enddo
           uc(im,j) = dtdxe(j)*(wk1(im,j)-wk1(1,j))  + dyce(j)*fy(im,j)+wku(im,j)
      enddo

#if defined(INNER_OMP)
!$omp parallel do default(shared) private(j,i)
#endif
      do j=js2g0,jn2g0
          do i=1,im
            vc(i,j) = dtdy*(wk1(i,j)-wk1(i,j+1)) - dx(j)*fx(i,j)+wkv(i,j)
          enddo
      enddo

 end subroutine d_sw
!----------------------------------------------------------------------- 

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: d2a2c_winds --- Interpolate winds
!
! !INTERFACE:
 subroutine d2a2c_winds(grid, u, v, ua, va, uc, vc, u_cen, v_cen, reset_winds, met_rlx)

  implicit none

! !PARAMETERS:
  type (T_FVDYCORE_GRID), intent(in) :: grid

  real(r8), intent(in   ):: u(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_s)
  real(r8), intent(inout):: v(grid%im,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d)
  real(r8), intent(  out):: ua(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
  real(r8), intent(  out):: va(grid%im,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d)
  real(r8), intent(  out):: uc(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
  real(r8), intent(  out):: vc(grid%im,grid%jfirst-2:grid%jlast+2 )

! cell centered winds
  logical , intent(in):: reset_winds
  real(r8), intent(in):: u_cen(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
  real(r8), intent(in):: v_cen(grid%im,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d)
  real(r8), intent(in):: met_rlx

! !DESCRIPTION:
!
!   Calculate the cell-centered (A-grid) winds and the cell-wall perpendicular
!   (C-grid) winds from the cell-wall parallel (D-grid) winds.
!
!   This routine assumes that U and V have complete haloes!  As a result,
!   the A-grid and C-grid results should have complete haloes from +/- ng_c
!   (which is generally smaller than ng_d).  This feature has not been 
!   thoroughly tested.
!
! !REVISION HISTORY:
!   WP   2007.06.01     Creation
!   WS   2007.07.03     Added ProTeX documentation, removed unused vars.
!   WS   2009.04.15     Fixed numerous problems in indexing bounds
!
!EOP
!-----------------------------------------------------------------------
!BOC

  real(r8) us, vs, un, vn
  real(r8) uanp(grid%im), uasp(grid%im), vanp(grid%im), vasp(grid%im), r2im

  real(r8), pointer:: sinlon(:)
  real(r8), pointer:: coslon(:)
  real(r8), pointer:: sinl5(:)
  real(r8), pointer:: cosl5(:)

  integer :: i, j, im2
  integer :: im, jm, jfirst, jlast, ng_s, ng_c, ng_d
  integer :: jn1gc, js1gc, jn2gc, js2gc   ! ng_c ghosted bounds
  integer :: js2gd, jn2gd                 ! ng_d ghosted bounds
  integer :: js2gs, jn2gsm1               ! ng_s ghosted bounds
  integer :: js2g2, jn1g2                 ! 2-lat ghosted bounds

  im     = grid%im
  jm     = grid%jm
  jfirst = grid%jfirst
  jlast  = grid%jlast

  ng_c   = grid%ng_c
  ng_d   = grid%ng_d
  ng_s   = grid%ng_s

  im2 = im/2

  js1gc  = max(1,jfirst-ng_c) ! ng_c lats on S (starting at 1)
  jn1gc  = min(jm,jlast+ng_c) ! ng_c latitudes on N (ending at jm)
  js2gc  = max(2,jfirst-ng_c) ! ng_c lats on S (starting at 2)
  jn2gc  = min(jm-1,jlast+ng_c)   ! ng_c latitudes on N (ending at jm-1)
  js2gs  = max(2,jfirst-ng_s) ! ng_s latitudes on S (starting at 2)
  jn2gsm1= min(jm-1,jlast+ng_s-1) ! ng_s-1 latitudes on N (ending at jm-1)
  js2gd  = max(2,jfirst-ng_d) ! ng_d latitudes on S (starting at 2)
  jn2gd  = min(jm-1,jlast+ng_d)   ! ng_d latitudes on N (ending at jm-1)
  js2g2  = max(2,jfirst-2)    ! 2 latitudes on S (starting at 2)
  jn1g2  = min(jm,jlast+2)    ! 2 latitudes on N (ending at jm)

  sinlon => grid%sinlon
  coslon => grid%coslon
  sinl5  => grid%sinl5
  cosl5  => grid%cosl5

! Get D-grid V-wind at the poles.

    r2im = 0.5_r16/real(im,r16)

    if ( jfirst-ng_d <= 1 ) then

!
! Treat SP
!
       do i=1,im-1
          uasp(i) = u(i,2) + u(i,3)
          vasp(i) = v(i,2) + v(i+1,2)
       enddo

       uasp(im) = u(im,2) + u(im,3)
       vasp(im) = v(im,2) + v(1,2)

! Projection at SP

       us = D0_0
       vs = D0_0

       do i=1,im2
          us = us + (uasp(i+im2)-uasp(i))*sinlon(i)    &
                  + (vasp(i)-vasp(i+im2))*coslon(i)
          vs = vs + (uasp(i+im2)-uasp(i))*coslon(i)    &
                  + (vasp(i+im2)-vasp(i))*sinlon(i)
       enddo
       us = us*r2im
       vs = vs*r2im

! get V-wind at SP

       do i=1,im2
          v(i,  1) =  us*cosl5(i) - vs*sinl5(i)
          v(i+im2,1) = -v(i,  1)
       enddo

    endif

    if ( jlast+ng_d >= jm ) then

!
! Treat NP
!
       do i=1,im-1
          uanp(i) = u(i,jm-1) + u(i,jm)
          vanp(i) = v(i,jm-1) + v(i+1,jm-1)
       enddo

       uanp(im) = u(im,jm-1) + u(im,jm)
       vanp(im) = v(im,jm-1) + v(1,jm-1)

! Projection at NP

       un = D0_0
       vn = D0_0
       do i=1,im2
          un = un + (uanp(i+im2)-uanp(i))*sinlon(i)   &
                  + (vanp(i+im2)-vanp(i))*coslon(i)
          vn = vn + (uanp(i)-uanp(i+im2))*coslon(i)   &
                  + (vanp(i+im2)-vanp(i))*sinlon(i)
       enddo
       un = un*r2im
       vn = vn*r2im

! get V-wind at NP

       do i=1,im2
          v(i,jm) = -un*cosl5(i) - vn*sinl5(i)
          v(i+im2,jm) = -v(i,jm)
       enddo

    endif

    ua(:,:) = D0_0
    va(:,:) = D0_0

    do j=js2gs,jn2gd
       do i=1,im-1
          va(i,j) = v(i,j) + v(i+1,j)
       enddo
          va(im,j) = v(im,j) + v(1,j)
    enddo

    do j=js2gd,jn2gsm1  ! WS: should be safe since u defined to jn2gs
       do i=1,im
          ua(i,j) = u(i,j) + u(i,j+1)
       enddo
    enddo

!
! reset cell center winds to the offline meteorlogy data
!

    if ( reset_winds .and. met_rlx > D0_0 ) then          
       ua(:,:) = (D1_0-met_rlx)*ua(:,:) + met_rlx*( D2_0*u_cen(:,:) )
       va(:,:) = (D1_0-met_rlx)*va(:,:) + met_rlx*( D2_0*v_cen(:,:) )
    endif

    if ( jfirst-ng_d <= 1 ) then
! Projection at SP
       us = D0_0
       vs = D0_0


       do i=1,im2
          us = us + (ua(i+im2,2)-ua(i    ,2))*sinlon(i)         &
               + (va(i    ,2)-va(i+im2,2))*coslon(i)
          vs = vs + (ua(i+im2,2)-ua(i    ,2))*coslon(i)         &
               + (va(i+im2,2)-va(i    ,2))*sinlon(i)
       enddo

       us = us/im
       vs = vs/im

       ! SP
       do i=1,im2
          ua(i,1)  = -us*sinlon(i) - vs*coslon(i)
          va(i,1)  =  us*coslon(i) - vs*sinlon(i)
          ua(i+im2,1)  = -ua(i,1)
          va(i+im2,1)  = -va(i,1)
       enddo

    endif

    if ( jlast+ng_d >= jm ) then

! Projection at NP
       un = D0_0
       vn = D0_0

       j = jm-1
       do i=1,im2
          un = un + (ua(i+im2,j)-ua(i    ,j))*sinlon(i)        &
               + (va(i+im2,j)-va(i    ,j))*coslon(i)
          vn = vn + (ua(i    ,j)-ua(i+im2,j))*coslon(i)        &
               + (va(i+im2,j)-va(i    ,j))*sinlon(i)
       enddo

       un = un/im
       vn = vn/im

       ! NP
       do i=1,im2
          ua(i,jm) = -un*sinlon(i) + vn*coslon(i)
          va(i,jm) = -un*coslon(i) - vn*sinlon(i)
          ua(i+im2,jm) = -ua(i,jm)
          va(i+im2,jm) = -va(i,jm)
       enddo

    endif


! A -> C
        do j=js1gc,jn1gc       ! uc needed N*ng_c S*ng_c, ua defined at poles
! i=1
            uc(1,j) = D0_25*(ua(1,j)+ua(im,j))


          do i=2,im

            uc(i,j) = D0_25*(ua(i,j)+ua(i-1,j))
          enddo
        enddo

        do j=js2g2,jn1g2       ! vc needed N*2, S*2 (for ycc), va defined at poles
          do i=1,im
            vc(i,j) = D0_25*(va(i,j)+va(i,j-1))  ! va needed N*2 S*3
          enddo
        enddo
        if ( jfirst==1 ) then
           do i=1,im
              vc(i,1) = D0_0   ! Needed to avoid undefined values in VC
           enddo
        endif
!EOC
 end subroutine d2a2c_winds
!-----------------------------------------------------------------------
 end module sw_core
