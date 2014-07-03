!=======================================================================
!BOP
!
! !MODULE: ice_transport_remap - horizontal transport via incremental remapping
!
! !DESCRIPTION:
!
! Transports quantities using the second-order conservative remapping
! scheme developed by John Dukowicz and John Baumgardner (DB) and modified
! for sea ice by William Lipscomb and Elizabeth Hunke.
!
! References:
!
! Dukowicz, J. K., and J. R. Baumgardner, 2000: Incremental
!  remapping as a transport/advection algorithm, J. Comput. Phys.,
!  160, 318-335.
!
! Lipscomb, W. H., and E. C. Hunke, 2004: Modeling sea ice
!  transport using incremental remapping, Mon. Wea. Rev., 132,
!  1341-1354.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_transport_remap.F 33 2006-11-13 19:51:14Z eclare $
!
! authors William H. Lipscomb, LANL
!         John Baumgardner, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004-05: Block structure added (WHL)
! 2006: Moved remap driver to ice_transport_driver
!       Geometry changes: 
!       (1) Reconstruct fields in stretched logically rectangular coordinates
!       (2) Modify geometry so that the area flux across each edge
!           can be specified (following an idea of Mats Bentsen)
!
! !INTERFACE:
!
      module ice_transport_remap
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
      use ice_domain_size
      use ice_constants
      use ice_fileunits, only: nu_diag, flush_fileunit
      use perf_mod,      only: t_startf, t_stopf, t_barrierf
!
!EOP
!
      implicit none
      save
      private
      public :: init_remap, horizontal_remap, make_masks

      logical (kind=log_kind), public :: maskhalo_remap

      integer (kind=int_kind), parameter ::                      &
         max_ntrace = 2+max_ntrcr+nilyr+nslyr  ! hice,hsno,qice,qsno,trcr

      integer (kind=int_kind), parameter ::     &
         ngroups  = 6      ,&! number of groups of triangles that
                             ! contribute transports across each edge
         nvert = 3           ! number of vertices in a triangle

      ! for triangle integral formulas
      real (kind=dbl_kind), parameter ::  & 
         p5625m = -9._dbl_kind/16._dbl_kind    ,&
         p52083 = 25._dbl_kind/48._dbl_kind

      logical (kind=log_kind), parameter :: bugcheck = .false.

!=======================================================================
! Here is some information about how the incremental remapping scheme
! works in CICE and how it can be adapted for use in other models.  
!
! The remapping routine is designed to transport a generic mass-like 
! field (in CICE, the ice fractional area) along with an arbitrary number
! of tracers in two dimensions.  The velocity components are assumed 
! to lie at grid cell corners and the transported scalars at cell centers. 
! Incremental remapping has the following desirable properties: 
! 
! (1) Tracer monotonicity is preserved.  That is, no new local 
!     extrema are produced in fields like ice thickness or internal 
!     energy. 
! (2) The reconstucted mass and tracer fields vary linearly in x and y. 
!     This means that remapping is 2nd-order accurate in space, 
!     except where horizontal gradients are limited to preserve 
!     monotonicity. 
! (3) There are economies of scale.  Transporting a single field 
!     is rather expensive, but additional fields have a relatively 
!     low marginal cost. 
! 
! The following generic conservation equations may be solved: 
! 
!            dm/dt = del*(u*m)             (0) 
!       d(m*T1)/dt = del*(u*m*T1)          (1) 
!    d(m*T1*T2)/dt = del*(u*m*T1*T2)       (2) 
! d(m*T1*T2*T3)/dt = del*(u*m*T1*T2*T3)    (3) 
!
! where d is a partial derivative, del is the 2D divergence operator,
! u is the horizontal velocity, m is the mass density field, and
! T1, T2, and T3 are tracers.
!
! In CICE, these equations have the form
! 
!               da/dt = del*(u*a)          (4)
! dv/dt =   d(a*h)/dt = del*(u*a*h)        (5)
! de/dt = d(a*h*q)/dt = del*(u*a*h*q)      (6)
!            d(aT)/dt = del*(u*a*t)        (7)
! 
! where a = fractional ice area, v = ice/snow volume, h = v/a = thickness, 
! e = ice/snow internal energy (J/m^2), q = e/v = internal energy per 
! unit volume (J/m^3), and T is a tracer.  These equations express 
! conservation of ice area, volume, internal energy, and area-weighted
! tracer, respectively. 
!
! (Note: In CICE, a, v and e are prognostic quantities from which
!  h and q are diagnosed.  The remapping routine works with tracers,
!  which means that h and q must be derived from a, v, and e before
!  calling the remapping routine.)  
!
! Earlier versions of CICE assumed fixed ice and snow density. 
! Beginning with CICE 4.0, the ice and snow density can be variable. 
! In this case, equations (5) and (6) are replaced by 
! 
! dv/dt =        d(a*h)/dt = del*(u*a*h)          (8)  
! dm/dt =    d(a*h*rho)/dt = del*(u*a*h*rho)      (9)
! de/dt = d(a*h*rho*qm)/dt = del*(u*a*h*rho*qm)   (10)
! 
! where rho = density and qm = internal energy per unit mass (J/kg). 
! Eq. (9) expresses mass conservation, which in the variable-density 
! case is no longer equivalent to volume conservation (8). 
!
! Tracers satisfying equations of the form (1) are called "type 1." 
! In CICE the paradigmatic type 1 tracers are hi and hs. 
! 
! Tracers satisfying equations of the form (2) are called "type 2". 
! The paradigmatic type 2 tracers are qi and qs (or rhoi and rhos 
!  in the variable-density case). 
! 
! Tracers satisfying equations of the form (3) are called "type 3."
! The paradigmatic type 3 tracers are qmi and qms in the variable-density
! case.  There are no such tracers in the constant-density case. 
! 
! The fields a, T1, and T2 are reconstructed in each grid cell with 
! 2nd-order accuracy.  T3 is reconstructed with 1st-order accuracy 
! (i.e., it is transported in upwind fashion) in order to avoid 
! additional mathematical complexity. 
! 
! The mass-like field lives in the array "mm" (shorthand for mean 
! mass) and the tracers fields in the array "tm" (mean tracers). 
! In order to transport tracers correctly, the remapping routine 
! needs to know the tracers types and relationships.  This is done 
! as follows: 
! 
! Each field in the "tm" array is assigned an index, 1:max_ntrace. 
! (Note: max_ntrace is not the same as max_ntrcr, the number of tracers 
! in the trcrn state variable array.  For remapping purposes we 
! have additional tracers hi, hs, qi and qs.) 
! For standard CICE with ntrcr = 1, nilyr = 4, and nslyr = 1, the 
! indexing is as follows: 
! 1   = hi 
! 2   = hs 
! 3   = Ts 
! 4-7 = qi 
! 8   = qs 
! 
! The tracer types (1,2,3) are contained in the "tracer_type" array. 
! For standard CICE: 
! 
!     tracer_type = (1 1 1 2 2 2 2 2) 
! 
! Type 2 and type 3 tracers are said to depend on type 1 tracers. 
! For instance, qi depends on hi, which is to say that 
! there is a conservation equation of the form (2) or (6). 
! Thus we define a "depend" array.  For standard CICE: 
! 
!          depend = (0 0 0 1 1 1 1 2) 
! 
! which implies that elements 1-3 (hi, hs, Ts) are type 1, 
! elements 4-7 (qi) depend on element 1 (hi), and element 8 (qs) 
! depends on element 2 (hs). 
!
! We also define a logical array "has_dependents".  In standard CICE: 
! 
!  has_dependents = (T T F F F F F F), 
! 
! which means that only elements 1 and 2 (hi and hs) have dependent 
! tracers. 
! 
! For the variable-density case, things are a bit more complicated. 
! Suppose we have 4 variable-density ice layers and one variable- 
! density snow layer.  Then the indexing is as follows: 
! 1    = hi 
! 2    = hs 
! 3    = Ts 
! 4-7  = rhoi 
! 8    = rhos 
! 9-12 = qmi 
! 13   = qms 
! 
! The key arrays are: 
! 
!    tracer_type = (1 1 1 2 2 2 2 2 3 3 3 3 3) 
! 
!         depend = (0 0 0 1 1 1 1 2 4 5 6 7 8) 
! 
! has_dependents = (T T F T T T T T F F F F F) 
! 
! which imply that hi and hs are type 1 with dependents rhoi and rhos, 
! while rhoi and rhos are type 2 with dependents qmi and qms. 
! 
! Tracers added to the ntrcr array are handled automatically 
! by the remapping with little extra coding.  It is necessary 
! only to provide the correct type and dependency information. 
!
! When using this routine in other models, most of the tracer dependency
! apparatus may be irrelevant.  In a layered ocean model, for example,
! the transported fields are the layer thickness h (the mass density
! field) and two or more tracers (T, S, and various trace species).
! Suppose there are just two tracers, T and S.  Then the tracer arrays
! have the values:
!
!    tracer_type = (1 1)
!         depend = (0 0)
! has_dependents = (F F)
!
! which is to say that all tracer transport equations are of the form (1).
!
! The tracer dependency arrays are optional input arguments for the
! main remapping subroutine.  If these arrays are not passed in, they
! take on the default values tracer_type(:) = 1, depend(:) = 0, and
! has_dependents(:) = F, which are appropriate for most purposes.
!
! Another optional argument is integral_order.  If integral_order = 1,
! then the triangle integrals are exact for linear functions of x and y.
! If integral_order = 2, these integrals are exact for both linear and
! quadratic functions.  If integral_order = 3, integrals are exact for
! cubic functions as well.  If all tracers are of type 1, then the
! integrals of mass*tracer are quadratic, and integral_order = 2 is
! sufficient.  In CICE, where there are type 2 tracers, we integrate
! functions of the form mass*tracer1*tracer2.  Thus integral_order = 3
! is required for exactness, though integral_order = 2 may be good enough
! in practice.
!
! Finally, a few words about the edgearea fields:
!
! In earlier versions of this scheme, the divergence of the velocity
! field implied by the remapping was, in general, different from the
! value of del*u computed in the dynamics.  For energetic consistency
! (in CICE as well as in layered ocean models such as HYPOP),
! these two values should agree.  This can be ensured by setting
! l_fixed_area = T and specifying the area transported across each grid
! cell edge in the arrays edgearea_e and edgearea_n.  The departure
! regions are then tweaked, following an idea by Mats Bentsen, such
! that they have the desired area.  If l_fixed_area = F, these regions
! are not tweaked, and the edgearea arrays are output variables.
!   
!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: init_remap - initialize grid quantities used for remapping
!
! !INTERFACE:
!
      subroutine init_remap
!
! !DESCRIPTION:
!
! Grid quantities used by the remapping transport scheme
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_blocks
      use ice_grid, only: dxt, dyt,                      &
                          xav, yav, xxav, xyav, yyav,    &
                          xxxav, xxyav, xyyav, yyyav
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::     &
 	 i, j, iblk     ! standard indices

      ! Compute grid cell average geometric quantities on the scaled
      ! rectangular grid with dx = 1, dy = 1.
      !
      ! Note: On a rectangular grid, the integral of any odd function
      !       of x or y = 0.

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            xav(i,j,iblk) = c0
            yav(i,j,iblk) = c0
!!!            These formulas would be used on a rectangular grid
!!!            with dimensions (dxt, dyt):  
!!!            xxav(i,j,iblk) = dxt(i,j,iblk)**2 / c12
!!!            yyav(i,j,iblk) = dyt(i,j,iblk)**2 / c12
            xxav(i,j,iblk) = c1/c12
            yyav(i,j,iblk) = c1/c12
            xyav(i,j,iblk) = c0
            xxxav(i,j,iblk) = c0
            xxyav(i,j,iblk) = c0
            xyyav(i,j,iblk) = c0
            yyyav(i,j,iblk) = c0
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
 
      end subroutine init_remap

!=======================================================================
!BOP
!
! !IROUTINE: horizontal_remap - incremental remapping transport scheme
!
! !INTERFACE:
!
      subroutine horizontal_remap (dt,                ntrace,     &
                                   uvel,              vvel,       &
                                   mm,                tm,         &
                                   l_fixed_area,                  &
                                   edgearea_e,        edgearea_n, &
                                   tracer_type_in,    depend_in,  &
                                   has_dependents_in,             &
                                   integral_order_in,             &
                                   l_dp_midpt_in)
!
! !DESCRIPTION:

! Solve the transport equations for one timestep using the incremental
! remapping scheme developed by John Dukowicz and John Baumgardner (DB)
! and modified for sea ice by William Lipscomb and Elizabeth Hunke.
!
! This scheme preserves monotonicity of ice area and tracers.  That is,
! it does not produce new extrema.  It is second-order accurate in space,
! except where gradients are limited to preserve monotonicity. 
!
! This version of the remapping allows the user to specify the areal
! flux across each edge, based on an idea developed by Mats Bentsen.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
! 2006: Moved driver (subroutine transport_remap) into separate module. 
!       Geometry changes (logically rectangular coordinates, fixed
!        area fluxes)
!       
! !USES:
!
      use ice_boundary
      use ice_global_reductions
      use ice_domain
      use ice_blocks
      use ice_grid, only: HTE, HTN, dxt, dyt, dxu, dyu,       &
                          tarea, tarear, hm,                  &
                          xav, yav, xxav, xyav, yyav,         &
                          xxxav, xxyav, xyyav, yyyav
      use ice_exit
      use ice_calendar, only: istep1
      use ice_timers
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::     &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         ntrace       ! number of tracers in use

      real (kind=dbl_kind), intent(in),       &
                dimension(nx_block,ny_block,max_blocks) ::           &
         uvel       ,&! x-component of velocity (m/s)
         vvel         ! y-component of velocity (m/s)

      real (kind=dbl_kind), intent(inout),     &
         dimension (nx_block,ny_block,0:ncat,max_blocks) ::          &
         mm           ! mean mass values in each grid cell

      real (kind=dbl_kind), intent(inout),     &
         dimension (nx_block,ny_block,max_ntrace,ncat,max_blocks) ::     &
         tm           ! mean tracer values in each grid cell

    !-------------------------------------------------------------------
    ! If l_fixed_area is true, the area of each departure region is
    !  computed in advance (e.g., by taking the divergence of the 
    !  velocity field and passed to locate_triangles.  The departure 
    !  regions are adjusted to obtain the desired area.
    ! If false, edgearea is computed in locate_triangles and passed out.
    !-------------------------------------------------------------------

      logical, intent(in) ::    &
         l_fixed_area     ! if true, edgearea_e and edgearea_n are prescribed
                          ! if false, edgearea is computed here and passed out

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks),  &
         intent(inout) ::                                             &
         edgearea_e     ,&! area of departure regions for east edges
         edgearea_n       ! area of departure regions for north edges

      integer (kind=int_kind), dimension (ntrace), intent(in),     &
         optional ::           &
         tracer_type_in       ,&! = 1, 2, or 3 (see comments above)
         depend_in              ! tracer dependencies (see above)

      logical (kind=log_kind), dimension (ntrace), intent(in),     &
         optional ::     &
         has_dependents_in      ! true if a tracer has dependent tracers

      integer (kind=int_kind), intent(in), optional ::     &
         integral_order_in      ! polynomial order for triangle integrals

      logical (kind=log_kind), intent(in), optional ::     &
         l_dp_midpt_in          ! if true, find departure points using
                                ! corrected midpoint velocity
!
!EOP
!
      ! local variables

      integer (kind=int_kind), dimension (ntrace) ::     &
         tracer_type       ,&! = 1, 2, or 3 (see comments above)
         depend              ! tracer dependencies (see above)

      logical (kind=log_kind), dimension (ntrace) ::     &
         has_dependents      ! true if a tracer has dependent tracers

      integer (kind=int_kind) ::     &
         integral_order      ! polynomial order for triangle integrals

      logical (kind=log_kind) ::     &
         l_dp_midpt          ! if true, find departure points using
                             ! corrected midpoint velocity

      integer (kind=int_kind) ::     &
         i, j           ,&! horizontal indices
         iblk           ,&! block indices
         ilo,ihi,jlo,jhi,&! beginning and end of physical domain
         n,              &! ice category index
         m                ! ice tracer index

      integer (kind=int_kind), dimension(0:ncat,max_blocks) ::     &
         icellsnc         ! number of cells with ice

      integer (kind=int_kind),     &
         dimension(nx_block*ny_block,0:ncat,max_blocks) ::     &
         indxinc, indxjnc   ! compressed i/j indices

      type (block) ::     &
         this_block       ! block information for current block

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) ::     &
         dpx            ,&! x coordinates of departure points at cell corners
         dpy              ! y coordinates of departure points at cell corners

      real (kind=dbl_kind), dimension(nx_block,ny_block,0:ncat,max_blocks) :: &
         mc             ,&! mass at geometric center of cell
         mx, my         ,&! limited derivative of mass wrt x and y
         mmask            ! = 1. if mass is present, = 0. otherwise

      real (kind=dbl_kind),      &
         dimension (nx_block,ny_block,max_ntrace,ncat,max_blocks) ::     &
         tc             ,&! tracer values at geometric center of cell
         tx, ty         ,&! limited derivative of tracer wrt x and y
         tmask            ! = 1. if tracer is present, = 0. otherwise

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat,max_blocks) ::     &
         mflxe, mflxn     ! mass transports across E and N cell edges

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrace,ncat,max_blocks) ::     &
         mtflxe, mtflxn   ! mass*tracer transports across E and N cell edges

      real (kind=dbl_kind), dimension (nx_block,ny_block,ngroups,max_blocks) ::     &
         triarea          ! area of east-edge departure triangle

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:nvert,ngroups,max_blocks) ::  &
         xp, yp           ! x and y coordinates of special triangle points
                          ! (need 4 points for triangle integrals)

      integer (kind=int_kind),     &
         dimension (nx_block,ny_block,ngroups,max_blocks) ::     &
         iflux          ,&! i index of cell contributing transport
         jflux            ! j index of cell contributing transport

      integer (kind=int_kind), dimension(ngroups,max_blocks) ::     &
         icellsng         ! number of cells with ice

      integer (kind=int_kind),     &
         dimension(nx_block*ny_block,ngroups,max_blocks) ::     &
         indxing, indxjng ! compressed i/j indices

      logical (kind=log_kind), dimension(max_blocks) ::     &
         l_stop           ! if true, abort the model

      integer (kind=int_kind), dimension(max_blocks) ::     &
         istop, jstop     ! indices of grid cell where model aborts

      character (len=char_len), dimension(max_blocks) ::   &
         edge             ! 'north' or 'east'

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka, &
         workb, &
         workc, &
         workd

      real (kind=dbl_kind), dimension (nx_block,ny_block,2,max_blocks) ::     &
          dpwork          
 
      real (kind=dbl_kind), dimension(nx_block,ny_block,2,0:ncat,max_blocks) :: &
          mwork

      integer (kind=int_kind), &
         dimension(nx_block,ny_block,max_blocks) :: halomask
      type (ice_halo) :: halo_info_tracer
    !------------------------------------------------------------------- 

      call t_barrierf('cice_hmap_remap1_BARRIER',MPI_COMM_ICE)
      call t_startf  ('cice_hmap_remap1')

      l_stop = .false.
      istop = 0
      jstop = 0

    !------------------------------------------------------------------- 
    ! Initialize various remapping arrays and options
    ! These are either passed in as optional arguments or set to the
    ! default values.
    !------------------------------------------------------------------- 

      if (present(tracer_type_in)) then
         tracer_type(:) = tracer_type_in(:)
      else
         tracer_type(:) = 1
      endif

      if (present(depend_in)) then
         depend(:) = depend_in(:)
      else
         depend(:) = 0
      endif

      if (present(has_dependents_in)) then
         has_dependents(:) = has_dependents_in(:)
      else
         has_dependents(:) = .false.
      endif

      if (present(integral_order_in)) then
         integral_order = integral_order_in
      else
         integral_order = 2   ! quadratic integrals
      endif

      if (present(l_dp_midpt_in)) then
         l_dp_midpt = l_dp_midpt_in
      else
         l_dp_midpt = .false.
      endif

      worka(:,:) = c1
      workb(:,:) = c1
      workc(:,:) = c1
      workd(:,:) = c1

!---!-------------------------------------------------------------------
!---! Remap the ice area and associated tracers.
!---! Remap the open water area (without tracers).
!---!-------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,n)
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

    !------------------------------------------------------------------- 
    ! Compute masks and count ice cells.
    ! Masks are used to prevent tracer values in cells without ice from
    !  being used to compute tracer gradients.
    !------------------------------------------------------------------- 

         call t_startf  ('cice_hmap_remap1_masks')
         call make_masks (nx_block,           ny_block,             &
                          ilo, ihi,           jlo, jhi,             &
                          nghost,             ntrace,               &
                          has_dependents,     icellsnc (:,iblk),    &
                          indxinc(:,:,iblk),  indxjnc  (:,:,iblk),  &
                          mm   (:,:,:,iblk),  mmask  (:,:,:,iblk),  &
                          tm (:,:,:,:,iblk),  tmask(:,:,:,:,iblk))
         call t_stopf  ('cice_hmap_remap1_masks')

    !-------------------------------------------------------------------
    ! Construct linear fields, limiting gradients to preserve monotonicity.
    ! Note: Pass in unit arrays instead of true distances HTE, HTN, etc.
    !       The resulting gradients are in scaled coordinates.
    !-------------------------------------------------------------------

         ! open water

         call t_startf  ('cice_hmap_remap1_cf1')
         call construct_fields(nx_block,           ny_block,           &
                               ilo, ihi,           jlo, jhi,           &
                               nghost,             ntrace,             &
                               tracer_type,        depend,             &
                               has_dependents,     icellsnc (0,iblk),  &
                               indxinc(:,0,iblk),  indxjnc(:,0,iblk),  &
!                               HTN    (:,:,iblk),  HTE    (:,:,iblk),  &
                               worka       (:,:),  workb       (:,:),  &
                               hm     (:,:,iblk),  xav    (:,:,iblk),  &
                               yav    (:,:,iblk),  xxav   (:,:,iblk),  &
                               xyav   (:,:,iblk),  yyav   (:,:,iblk),  &
                               xxxav  (:,:,iblk),  xxyav  (:,:,iblk),  &
                               xyyav  (:,:,iblk),  yyyav  (:,:,iblk),  &
!                               dxt    (:,:,iblk),  dyt    (:,:,iblk),  &
                               workc       (:,:),  workd       (:,:),  &
                               mm   (:,:,0,iblk),  mc   (:,:,0,iblk),  &
                               mx   (:,:,0,iblk),  my   (:,:,0,iblk),  &
                               mmask(:,:,0,iblk) )
         call t_stopf  ('cice_hmap_remap1_cf1')

         ! ice categories

         call t_startf  ('cice_hmap_remap1_cf2')
         do n = 1, ncat

            call construct_fields(nx_block,             ny_block,           &
                                  ilo, ihi,             jlo, jhi,           &
                                  nghost,               ntrace,             &
                                  tracer_type,          depend,             &
                                  has_dependents,       icellsnc (n,iblk),  &
                                  indxinc  (:,n,iblk),  indxjnc(:,n,iblk),  &
!                                  HTN      (:,:,iblk),  HTE    (:,:,iblk),  &
                                  worka         (:,:),  workb       (:,:),  &
                                  hm       (:,:,iblk),  xav    (:,:,iblk),  &
                                  yav      (:,:,iblk),  xxav   (:,:,iblk),  &
                                  xyav     (:,:,iblk),  yyav   (:,:,iblk),  &
                                  xxxav    (:,:,iblk),  xxyav  (:,:,iblk),  &
                                  xyyav    (:,:,iblk),  yyyav  (:,:,iblk),  &
!                                  dxt      (:,:,iblk),  dyt   (:,:,iblk),   &
                                  workc         (:,:),  workd       (:,:),  &
                                  mm     (:,:,n,iblk),  mc   (:,:,n,iblk),  &
                                  mx     (:,:,n,iblk),  my   (:,:,n,iblk),  &
                                  mmask  (:,:,n,iblk),                      &
                                  tm   (:,:,:,n,iblk),  tc (:,:,:,n,iblk),  &
                                  tx   (:,:,:,n,iblk),  ty (:,:,:,n,iblk),  &
                                  tmask(:,:,:,n,iblk) )

         enddo                  ! n
         call t_stopf  ('cice_hmap_remap1_cf2')

    !-------------------------------------------------------------------
    ! Given velocity field at cell corners, compute departure points
    ! of trajectories.
    !-------------------------------------------------------------------

         call t_startf  ('cice_hmap_remap1_dp')
         call departure_points(nx_block,        ny_block,        &
                               ilo, ihi,        jlo, jhi,        &
                               nghost,          dt,              &
                               uvel(:,:,iblk),  vvel(:,:,iblk),  &
                               dxu (:,:,iblk),  dyu (:,:,iblk),  &
                               HTN (:,:,iblk),  HTE (:,:,iblk),  &
                               dpx (:,:,iblk),  dpy (:,:,iblk),  &
                               l_dp_midpt,      l_stop  (iblk),  &
                               istop   (iblk),  jstop   (iblk))
         call t_stopf  ('cice_hmap_remap1_dp')

         if (l_stop(iblk)) then
            this_block = get_block(blocks_ice(iblk),iblk)         
            write(nu_diag,*) 'istep1, my_task, iblk =',            &
                              istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop(iblk) > 0 .and. jstop(iblk) > 0)             &
                 write(nu_diag,*) 'Global i and j:',               &
                                  this_block%i_glob(istop(iblk)),  &
                                  this_block%j_glob(jstop(iblk)) 
            call abort_ice('remap transport: bad departure points')
         endif

      enddo                     ! iblk
      !$OMP END PARALLEL DO

     call t_stopf ('cice_hmap_remap1')
    !-------------------------------------------------------------------
    ! Ghost cell updates
    ! If nghost >= 2, these calls are not needed
    !-------------------------------------------------------------------

      if (nghost==1) then

         call t_barrierf ('cice_hmap_halo1_BARRIER',MPI_COMM_ICE)
         call t_startf ('cice_hmap_halo1')
         call ice_timer_start(timer_bound)

         ! departure points

         ! load em up
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
              dpwork(i,j,1,  iblk) = dpx(i,j,iblk)
              dpwork(i,j,2,  iblk) = dpy(i,j,iblk)
              mwork (i,j,1,:,iblk) = mx (i,j,:,iblk)
              mwork (i,j,2,:,iblk) = my (i,j,:,iblk)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
!jw         call ice_HaloUpdate (dpx,                halo_info, &
!jw                              field_loc_NEcorner, field_type_vector)
!jw         call ice_HaloUpdate (dpy,                halo_info, &
!jw                              field_loc_NEcorner, field_type_vector)
         call ice_HaloUpdate (dpwork,             halo_info, &
                              field_loc_NEcorner, field_type_vector)

         ! mass field
         call ice_HaloUpdate (mc,               halo_info, &
                              field_loc_center, field_type_scalar)
!jw         call ice_HaloUpdate (mx,               halo_info, &
!jw                              field_loc_center, field_type_vector)
!jw         call ice_HaloUpdate (my,               halo_info, &
!jw                              field_loc_center, field_type_vector)
         call ice_HaloUpdate (mwork,            halo_info, &
                              field_loc_center, field_type_vector)

        call t_stopf ('cice_hmap_halo1')
        call t_barrierf ('cice_hmap_copy1_BARRIER',MPI_COMM_ICE)
        call t_startf ('cice_hmap_copy1')

         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
              dpx(i,j,  iblk) = dpwork(i,j,1,  iblk)
              dpy(i,j,  iblk) = dpwork(i,j,2,  iblk)
              mx (i,j,:,iblk) = mwork (i,j,1,:,iblk)
              my (i,j,:,iblk) = mwork (i,j,2,:,iblk)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         call t_stopf ('cice_hmap_copy1')

         call t_barrierf ('cice_hmap_halo2_BARRIER',MPI_COMM_ICE)
         call t_startf ('cice_hmap_halo2')

         if (maskhalo_remap) then
            call t_startf ('cice_hmap_halo2hm')
            halomask = 0
            !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,n,m,j,i)
            do iblk = 1,nblocks
               this_block = get_block(blocks_ice(iblk),iblk)         
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do n = 1,ncat
               do m = 1,ntrace
               do j = jlo, jhi
               do i = ilo, ihi
                  if (tc(i,j,m,n,iblk) /= c0) halomask(i,j,iblk) = 1
                  if (tx(i,j,m,n,iblk) /= c0) halomask(i,j,iblk) = 1
                  if (ty(i,j,m,n,iblk) /= c0) halomask(i,j,iblk) = 1
               enddo
               enddo
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
            call t_stopf ('cice_hmap_halo2hm')

            call t_barrierf ('cice_hmap_halo2hh_BARRIER',MPI_COMM_ICE)
            call t_startf ('cice_hmap_halo2hh')
            call ice_HaloUpdate(halomask, halo_info, &
                                field_loc_center, field_type_scalar)
            call t_stopf ('cice_hmap_halo2hh')

            call t_barrierf ('cice_hmap_halo2hc_BARRIER',MPI_COMM_ICE)
            call t_startf ('cice_hmap_halo2hc')
            call ice_HaloMask(halo_info_tracer, halo_info, halomask)
            call t_stopf ('cice_hmap_halo2hc')

            call t_barrierf ('cice_hmap_halo2hu_BARRIER',MPI_COMM_ICE)
            call t_startf ('cice_hmap_halo2hu')

            ! tracer fields 
            call ice_HaloUpdate (tc(:,:,1:ntrace,:,:), halo_info_tracer, &
                                 field_loc_center, field_type_scalar)
            call ice_HaloUpdate (tx(:,:,1:ntrace,:,:), halo_info_tracer, &
                                 field_loc_center, field_type_vector)
            call ice_HaloUpdate (ty(:,:,1:ntrace,:,:), halo_info_tracer, &
                                 field_loc_center, field_type_vector)
            call ice_timer_stop(timer_bound)
            call t_stopf ('cice_hmap_halo2hu')

            call t_barrierf ('cice_hmap_halo2hd_BARRIER',MPI_COMM_ICE)
            call t_startf ('cice_hmap_halo2hd')
            call ice_HaloDestroy(halo_info_tracer)
            call t_stopf ('cice_hmap_halo2hd')
         else
            ! tracer fields 
            call ice_HaloUpdate (tc(:,:,1:ntrace,:,:), halo_info, &
                                 field_loc_center, field_type_scalar)
            call ice_HaloUpdate (tx(:,:,1:ntrace,:,:), halo_info, &
                                 field_loc_center, field_type_vector)
            call ice_HaloUpdate (ty(:,:,1:ntrace,:,:), halo_info, &
                                 field_loc_center, field_type_vector)
         endif

         call ice_timer_stop(timer_bound)
         call t_stopf ('cice_hmap_halo2')

      endif  ! nghost

      call t_barrierf ('cice_hmap_remap2_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_hmap_remap2')

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,n)
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

    !-------------------------------------------------------------------
    ! Transports for east cell edges.
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Compute areas and vertices of departure triangles.
    !-------------------------------------------------------------------

         call t_startf  ('cice_hmap_remap2e')
         edge(iblk) = 'east'
         call locate_triangles(nx_block,              ny_block,           &
                               ilo, ihi,              jlo, jhi,           &
                               nghost,                edge(iblk),         &
                               icellsng    (:,iblk),                      &
                               indxing   (:,:,iblk),  indxjng(:,:,iblk),  &
                               dpx       (:,:,iblk),  dpy    (:,:,iblk),  &
                               dxu       (:,:,iblk),  dyu    (:,:,iblk),  &
                               xp    (:,:,:,:,iblk),  yp (:,:,:,:,iblk),  &
                               iflux   (:,:,:,iblk),  jflux(:,:,:,iblk),  &
                               triarea (:,:,:,iblk),  l_fixed_area,       &
                               edgearea_e(:,:,iblk))

    !-------------------------------------------------------------------
    ! Given triangle vertices, compute coordinates of triangle points
    !  needed for transport integrals.
    !-------------------------------------------------------------------

         call triangle_coordinates (nx_block,           ny_block,           &
                                    integral_order,     icellsng (:,iblk),  &
                                    indxing(:,:,iblk),  indxjng(:,:,iblk),  &
                                    xp (:,:,:,:,iblk),  yp (:,:,:,:,iblk))

    !-------------------------------------------------------------------
    ! Compute the transport across east cell edges by summing contributions
    ! from each triangle.
    !-------------------------------------------------------------------

         ! open water

         call transport_integrals(nx_block,           ny_block,             &
                                  ntrace,             icellsng (:,iblk),    &
                                  indxing(:,:,iblk),  indxjng  (:,:,iblk),  &
                                  tracer_type,        depend,               &
                                  integral_order,     triarea(:,:,:,iblk),  &
                                  iflux(:,:,:,iblk),  jflux  (:,:,:,iblk),  &
                                  xp (:,:,:,:,iblk),  yp   (:,:,:,:,iblk),  &
                                  mc   (:,:,0,iblk),  mx     (:,:,0,iblk),  &
                                  my   (:,:,0,iblk),  mflxe  (:,:,0,iblk))

         ! ice categories
         do n = 1, ncat
            call transport_integrals                                       &
                               (nx_block,           ny_block,              &
                                ntrace,             icellsng (:,iblk),     &
                                indxing(:,:,iblk),  indxjng   (:,:,iblk),  &
                                tracer_type,        depend,                &
                                integral_order,     triarea (:,:,:,iblk),  &
                                iflux(:,:,:,iblk),  jflux   (:,:,:,iblk),  &
                                xp (:,:,:,:,iblk),  yp    (:,:,:,:,iblk),  &
                                mc   (:,:,n,iblk),  mx      (:,:,n,iblk),  &
                                my   (:,:,n,iblk),  mflxe   (:,:,n,iblk),  &
                                tc (:,:,:,n,iblk),  tx    (:,:,:,n,iblk),  &
                                ty (:,:,:,n,iblk),  mtflxe(:,:,:,n,iblk))

         enddo
         call t_stopf  ('cice_hmap_remap2e')

    !-------------------------------------------------------------------
    ! Repeat for north edges
    !-------------------------------------------------------------------

         call t_startf  ('cice_hmap_remap2n')
         edge(iblk) = 'north'
         call locate_triangles(nx_block,              ny_block,           &
                               ilo, ihi,              jlo, jhi,           &
                               nghost,                edge(iblk),         &
                               icellsng    (:,iblk),                      &
                               indxing   (:,:,iblk),  indxjng(:,:,iblk),  &
                               dpx       (:,:,iblk),  dpy    (:,:,iblk),  &
                               dxu       (:,:,iblk),  dyu    (:,:,iblk),  &
                               xp    (:,:,:,:,iblk),  yp (:,:,:,:,iblk),  &
                               iflux   (:,:,:,iblk),  jflux(:,:,:,iblk),  &
                               triarea (:,:,:,iblk),  l_fixed_area,       &
                               edgearea_n(:,:,iblk))

         call triangle_coordinates (nx_block,           ny_block,           &
                                    integral_order,     icellsng (:,iblk),  &
                                    indxing(:,:,iblk),  indxjng(:,:,iblk),  &
                                    xp (:,:,:,:,iblk),  yp (:,:,:,:,iblk))

         ! open water
         call transport_integrals(nx_block,           ny_block,             &
                                  ntrace,             icellsng (:,iblk),    &
                                  indxing(:,:,iblk),  indxjng(:,:,iblk),    &
                                  tracer_type,        depend,               &
                                  integral_order,     triarea(:,:,:,iblk),  &
                                  iflux(:,:,:,iblk),  jflux  (:,:,:,iblk),  &
                                  xp (:,:,:,:,iblk),  yp   (:,:,:,:,iblk),  &
                                  mc   (:,:,0,iblk),  mx     (:,:,0,iblk),  &
                                  my   (:,:,0,iblk),  mflxn  (:,:,0,iblk))

         ! ice categories
         do n = 1, ncat
            call transport_integrals                                       &
                               (nx_block,           ny_block,              &
                                ntrace,             icellsng (:,iblk),     &
                                indxing(:,:,iblk),  indxjng   (:,:,iblk),  &
                                tracer_type,        depend,                &
                                integral_order,     triarea (:,:,:,iblk),  &
                                iflux(:,:,:,iblk),  jflux   (:,:,:,iblk),  &
                                xp (:,:,:,:,iblk),  yp    (:,:,:,:,iblk),  &
                                mc   (:,:,n,iblk),  mx      (:,:,n,iblk),  &
                                my   (:,:,n,iblk),  mflxn   (:,:,n,iblk),  &
                                tc (:,:,:,n,iblk),  tx    (:,:,:,n,iblk),  &
                                ty (:,:,:,n,iblk),  mtflxn(:,:,:,n,iblk))

         enddo                  ! n
         call t_stopf  ('cice_hmap_remap2n')

    !-------------------------------------------------------------------
    ! Update the ice area and tracers.
    !-------------------------------------------------------------------

         ! open water

         call t_startf  ('cice_hmap_remap2_upd')
         call update_fields (nx_block,           ny_block,          &
                             ilo, ihi,           jlo, jhi,          &
                             ntrace,                                &
                             tracer_type,        depend,            &
                             tarear(:,:,iblk),   l_stop(iblk),      &
                             istop(iblk),        jstop (iblk),      &
                             mflxe(:,:,0,iblk),  mflxn(:,:,0,iblk), &
                             mm   (:,:,0,iblk))

         if (l_stop(iblk)) then
            this_block = get_block(blocks_ice(iblk),iblk)         
            write (nu_diag,*) 'istep1, my_task, iblk, cat =',      &
                               istep1, my_task, iblk, '0'
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop(iblk) > 0 .and. jstop(iblk) > 0)             &
                 write(nu_diag,*) 'Global i and j:',               &
                                  this_block%i_glob(istop(iblk)),  &
                                  this_block%j_glob(jstop(iblk)) 
            call abort_ice ('ice remap_transport: negative area (open water)')
         endif


         ! ice categories
         do n = 1, ncat

            call update_fields(nx_block,              ny_block,              &
                               ilo, ihi,              jlo, jhi,              &
                               ntrace,                                       &
                               tracer_type,           depend,                &
                               tarear(:,:,iblk),      l_stop(iblk),          &
                               istop(iblk),           jstop(iblk),           &
                               mflxe (:,:,  n,iblk),  mflxn (:,:,  n,iblk),  &
                               mm    (:,:,  n,iblk),                         &
                               mtflxe(:,:,:,n,iblk),  mtflxn(:,:,:,n,iblk),  &
                               tm    (:,:,:,n,iblk))

            if (l_stop(iblk)) then
               this_block = get_block(blocks_ice(iblk),iblk)         
               write (nu_diag,*) 'istep1, my_task, iblk, cat =',      &
                                  istep1, my_task, iblk, n
               write (nu_diag,*) 'Global block:', this_block%block_id
               if (istop(iblk) > 0 .and. jstop(iblk) > 0)             &
                    write(nu_diag,*) 'Global i and j:',               &
                                     this_block%i_glob(istop(iblk)),  &
                                     this_block%j_glob(jstop(iblk)) 
               call abort_ice ('ice remap_transport: negative area (ice)')
            endif
         enddo                  ! n
         call t_stopf  ('cice_hmap_remap2_upd')

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call t_stopf ('cice_hmap_remap2')

      end subroutine horizontal_remap

!=======================================================================
!
!BOP
!
! !IROUTINE: make_masks - make area and tracer masks
!
! !INTERFACE:
!
      subroutine make_masks (nx_block, ny_block,           &
                             ilo, ihi, jlo, jhi,           &
                             nghost,   ntrace,             &
                             has_dependents,               &
                             icells,                       &
                             indxi,    indxj,              &
                             mm,       mmask,              &
                             tm,       tmask)

!
! !DESCRIPTION:
!
! Make area and tracer masks.
!
! If an area is masked out (mm < puny), then the values of tracers
!  in that grid cell are assumed to have no physical meaning.
!
! Similarly, if a tracer with dependents is masked out
!  (abs(tm) < puny), then the values of its dependent tracers in that
!  grid cell are assumed to have no physical meaning.
! For example, the enthalpy value has no meaning if the thickness
!  is zero.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block  ,&! block dimensions
           ilo,ihi,jlo,jhi     ,&! beginning and end of physical domain
           nghost              ,&! number of ghost cells
           ntrace                ! number of tracers in use

      logical (kind=log_kind), dimension (ntrace), intent(in) ::     &
           has_dependents      ! true if a tracer has dependent tracers

      integer (kind=int_kind), dimension(0:ncat), intent(out) ::     &
           icells         ! number of cells with ice

      integer (kind=int_kind), dimension(nx_block*ny_block,0:ncat),     &
           intent(out) ::     &
           indxi        ,&! compressed i/j indices
           indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
           intent(in) ::     &
           mm            ! mean ice area in each grid cell

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
           intent(out) ::     &
           mmask         ! = 1. if ice is present, else = 0.

      real (kind=dbl_kind), dimension (nx_block, ny_block, max_ntrace, ncat),  &
           intent(in), optional ::     &
           tm            ! mean tracer values in each grid cell

      real (kind=dbl_kind), dimension (nx_block, ny_block, max_ntrace, ncat),  &
           intent(out), optional ::     &
           tmask         ! = 1. if tracer is present, else = 0.
!
!EOP
!
      integer (kind=int_kind) ::     &
           i, j, ij       ,&! horizontal indices
           n              ,&! ice category index
           nt               ! tracer index

      do n = 0, ncat
         do ij = 1, nx_block*ny_block
            indxi(ij,n) = 0
            indxj(ij,n) = 0
         enddo
      enddo

    !-------------------------------------------------------------------
    ! open water mask
    !-------------------------------------------------------------------

      icells(0) = 0
      do j = 1, ny_block
      do i = 1, nx_block
         if (mm(i,j,0) > puny) then
            mmask(i,j,0) = c1
            icells(0) = icells(0) + 1
            ij = icells(0)
            indxi(ij,0) = i
            indxj(ij,0) = j
         else
            mmask(i,j,0) = c0
         endif
      enddo
      enddo

      do n = 1, ncat

    !-------------------------------------------------------------------
    ! Find grid cells where ice is present.
    !-------------------------------------------------------------------

         icells(n) = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (mm(i,j,n) > puny) then
               icells(n) = icells(n) + 1
               ij = icells(n)
               indxi(ij,n) = i
               indxj(ij,n) = j
            endif               ! mm > puny
         enddo
         enddo

    !-------------------------------------------------------------------
    ! ice area mask
    !-------------------------------------------------------------------

         mmask(:,:,n) = c0
         do ij = 1, icells(n)
            i = indxi(ij,n)
            j = indxj(ij,n)
            mmask(i,j,n) = c1
         enddo

    !-------------------------------------------------------------------
    ! tracer masks
    !-------------------------------------------------------------------

         if (present(tm)) then

            tmask(:,:,:,n) = c0
            do nt = 1, ntrace
               if (has_dependents(nt)) then
                  do ij = 1, icells(n)
                     i = indxi(ij,n)
                     j = indxj(ij,n)
                     if (abs(tm(i,j,nt,n)) > puny) then
                        tmask(i,j,nt,n) = c1
                     endif
                  enddo
               endif
            enddo

         endif                     ! present(tm)

    !-------------------------------------------------------------------
    ! Redefine icells
    ! For nghost = 1, exclude ghost cells
    ! For nghost = 2, include one layer of ghost cells
    !-------------------------------------------------------------------

         icells(n) = 0
         do j = jlo-nghost+1, jhi+nghost-1
         do i = ilo-nghost+1, ihi+nghost-1
            if (mm(i,j,n) > puny) then
               icells(n) = icells(n) + 1
               ij = icells(n)
               indxi(ij,n) = i
               indxj(ij,n) = j
            endif               ! mm > puny
         enddo
         enddo
      
      enddo ! n

      end subroutine make_masks

!=======================================================================
!
!BOP
!
! !IROUTINE: construct_fields - construct fields of ice area and tracers
!
! !INTERFACE:
!
      subroutine construct_fields (nx_block,       ny_block,   &
                                   ilo, ihi,       jlo, jhi,   &
                                   nghost,         ntrace,     &
                                   tracer_type,    depend,     &
                                   has_dependents, icells,     &
                                   indxi,          indxj,      &
                                   HTN,            HTE,        &
                                   hm,             xav,        &
                                   yav,            xxav,       &
                                   xyav,           yyav,       &
                                   xxxav,          xxyav,      &
                                   xyyav,          yyyav,      &
                                   dxt,            dyt,        &
                                   mm,             mc,         &
                                   mx,             my,         &
                                   mmask,                      &
                                   tm,             tc,         &
                                   tx,             ty,         &
                                   tmask)
!
! !DESCRIPTION:
!
! Construct fields of ice area and tracers.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         nx_block, ny_block  ,&! block dimensions
         ilo,ihi,jlo,jhi     ,&! beginning and end of physical domain
         nghost              ,&! number of ghost cells
         ntrace              ,&! number of tracers in use
         icells                ! number of cells with mass

      integer (kind=int_kind), dimension (ntrace), intent(in) ::     &
         tracer_type       ,&! = 1, 2, or 3 (see comments above)
         depend              ! tracer dependencies (see above)

      logical (kind=log_kind), dimension (ntrace), intent(in) ::     &
         has_dependents      ! true if a tracer has dependent tracers

      integer (kind=int_kind), dimension(nx_block*ny_block), intent(in) :: &
         indxi          ,&! compressed i/j indices
         indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block),   &
         intent(in) ::   &
         hm             ,&! land/boundary mask, thickness (T-cell)
         HTN            ,&! length of northern edge of T-cell (m)
         HTE            ,&! length of eastern edge of T-cell (m)
         xav,  yav              ,&! mean T-cell values of x, y
         xxav, xyav, yyav       ,&! mean T-cell values of xx, xy, yy
         xxxav,xxyav,xyyav,yyyav,&! mean T-cell values of , xxy, xyy, yyy
         dxt            ,&! grid cell width (m)
         dyt              ! grid cell height (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block),   &
         intent(in) ::   &
         mm            ,&! mean value of mass field
         mmask           ! = 1. if ice is present, = 0. otherwise

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrace),   &
         intent(in), optional ::   &
         tm             ,&! mean tracer
         tmask            ! = 1. if tracer is present, = 0. otherwise

      real (kind=dbl_kind), dimension (nx_block,ny_block),   &
         intent(out) ::   &
         mc             ,&! mass value at geometric center of cell
         mx, my           ! limited derivative of mass wrt x and y

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrace),   &
         intent(out), optional ::   &
         tc             ,&! tracer at geometric center of cell
         tx, ty           ! limited derivative of tracer wrt x and y
!
!EOP
!
      integer (kind=int_kind) ::   &
         i, j           ,&! horizontal indices
         nt, nt1        ,&! tracer indices
         ij               ! combined i/j horizontal index

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::    &
         mxav           ,&! x coordinate of center of mass
         myav             ! y coordinate of center of mass

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrace) ::  &
         mtxav          ,&! x coordinate of center of mass*tracer
         mtyav            ! y coordinate of center of mass*tracer

      real (kind=dbl_kind) ::   &
         w1, w2, w3, w4, w5, w6, w7   ! work variables

    !-------------------------------------------------------------------
    ! Compute field values at the geometric center of each grid cell,
    ! and compute limited gradients in the x and y directions.
    !
    ! For second order accuracy, each state variable is approximated as
    ! a field varying linearly over x and y within each cell.  For each
    ! category, the integrated value of m(x,y) over the cell must
    ! equal mm(i,j,n)*tarea(i,j), where tarea(i,j) is the cell area.
    ! Similarly, the integrated value of m(x,y)*t(x,y) must equal
    ! the total mass*tracer, mm(i,j,n)*tm(i,j,n)*tarea(i,j).
    !
    ! These integral conditions are satisfied for linear fields if we
    ! stipulate the following:
    ! (1) The mean mass, mm, is equal to the mass at the cell centroid.
    ! (2) The mean value tm1 of type 1 tracers is equal to the value
    !     at the center of mass.
    ! (3) The mean value tm2 of type 2 tracers is equal to the value
    !     at the center of mass*tm1, where tm2 depends on tm1.
    !     (See comments at the top of the module.)
    !
    ! We want to find the value of each state variable at a standard
    ! reference point, which we choose to be the geometric center of
    ! the cell.  The geometric center is located at the intersection
    ! of the line joining the midpoints of the north and south edges
    ! with the line joining the midpoints of the east and west edges.
    ! To find the value at the geometric center, we must know the
    ! location of the cell centroid/center of mass, along with the
    ! mean value and the gradients with respect to x and y.
    !
    ! The cell gradients are first computed from the difference between
    ! values in the neighboring cells, then limited by requiring that
    ! no new extrema are created within the cell.
    !
    ! For rectangular coordinates the centroid and the geometric
    ! center coincide, which means that some of the equations in this
    ! subroutine could be simplified.  However, the full equations
    ! are retained for generality.
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         mc(i,j)  = c0
         mx(i,j)  = c0
         my(i,j)  = c0
         mxav(i,j) = c0
         myav(i,j) = c0
      enddo
      enddo

      if (present(tm)) then
         do nt = 1, ntrace
            do j = 1, ny_block
            do i = 1, nx_block
               tc(i,j,nt) = c0
               tx(i,j,nt) = c0
               ty(i,j,nt) = c0
            enddo
            enddo
         enddo
      endif
         
      ! limited gradient of mass field in each cell (except masked cells)
      ! Note: The gradient is computed in scaled coordinates with
      !       dxt = dyt = hte = htn = 1.

      call limited_gradient (nx_block, ny_block,   &
                             ilo, ihi, jlo, jhi,   &
                             nghost,               &
                             mm,       hm,         &
                             xav,      yav,        &
                             HTN,      HTE,        &
                             dxt,      dyt,        &
                             mx,       my)

      do ij = 1,icells   ! ice is present
         i = indxi(ij)
         j = indxj(ij)

         ! mass field at geometric center
         mc(i,j) = mm(i,j) - xav(i,j)*mx(i,j)   &
                           - yav(i,j)*my(i,j)

      enddo                     ! ij

      ! tracers

      if (present(tm)) then

       do ij = 1,icells       ! cells with mass
          i = indxi(ij)
          j = indxj(ij)

         ! center of mass (mxav,myav) for each cell
          mxav(i,j) = (mx(i,j)*xxav(i,j)    &
                     + my(i,j)*xyav(i,j)    &
                     + mc(i,j)*xav (i,j)) / mm(i,j)
          myav(i,j) = (mx(i,j)*xyav(i,j)    &
                     + my(i,j)*yyav(i,j)    &
                     + mc(i,j)*yav(i,j)) / mm(i,j)
       enddo

       do nt = 1, ntrace

         if (tracer_type(nt)==1) then ! independent of other tracers

            call limited_gradient(nx_block,     ny_block,  &
                                  ilo, ihi,     jlo, jhi,  &
                                  nghost,                  &
                                  tm(:,:,nt),   mmask,     &
                                  mxav,         myav,      &
                                  HTN,          HTE,       &
                                  dxt,          dyt,       &
                                  tx(:,:,nt),   ty(:,:,nt)) 

            if (has_dependents(nt)) then   ! need center of area*tracer

               do j = 1, ny_block
               do i = 1, nx_block
                  mtxav(i,j,nt) = c0
                  mtyav(i,j,nt) = c0
               enddo
               enddo

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells  ! Note: no tx or ty in ghost cells
                                  ! (bound calls are later)
                  i = indxi(ij)
                  j = indxj(ij)

                  ! tracer value at geometric center
                  tc(i,j,nt) = tm(i,j,nt) - tx(i,j,nt)*mxav(i,j)   &
                                          - ty(i,j,nt)*myav(i,j)

                  if (tmask(i,j,nt) > puny) then

                     ! center of area*tracer
                     w1 = mc(i,j)*tc(i,j,nt)
                     w2 = mc(i,j)*tx(i,j,nt)   &
                        + mx(i,j)*tc(i,j,nt)
                     w3 = mc(i,j)*ty(i,j,nt)   &
                        + my(i,j)*tc(i,j,nt)
                     w4 = mx(i,j)*tx(i,j,nt)
                     w5 = mx(i,j)*ty(i,j,nt)   &
                        + my(i,j)*tx(i,j,nt)
                     w6 = my(i,j)*ty(i,j,nt)
                     w7 = c1 / (mm(i,j)*tm(i,j,nt))
                     mtxav(i,j,nt) = (w1*xav (i,j)  + w2*xxav (i,j)   &
                                    + w3*xyav (i,j) + w4*xxxav(i,j)   &
                                    + w5*xxyav(i,j) + w6*xyyav(i,j))   &
                                    * w7
                     mtyav(i,j,nt) = (w1*yav(i,j)   + w2*xyav (i,j)   &
                                    + w3*yyav(i,j)  + w4*xxyav(i,j)   &
                                    + w5*xyyav(i,j) + w6*yyyav(i,j))   &
                                    * w7
                  endif         ! tmask

               enddo            ! ij

            else                ! no dependents

               do ij = 1, icells      ! mass is present
                  i = indxi(ij)
                  j = indxj(ij)

                  ! tracer value at geometric center
                  tc(i,j,nt) = tm(i,j,nt) - tx(i,j,nt)*mxav(i,j)   &
                                          - ty(i,j,nt)*myav(i,j)
               enddo            ! ij

            endif               ! has_dependents

         elseif (tracer_type(nt)==2) then   ! tracer nt depends on nt1
            nt1 = depend(nt)

            call limited_gradient(nx_block,       ny_block,         &
                                  ilo, ihi,       jlo, jhi,         &
                                  nghost,                           &
                                  tm(:,:,nt),     tmask(:,:,nt1),   &
                                  mtxav(:,:,nt1), mtyav(:,:,nt1),   &
                                  HTN,            HTE,              &
                                  dxt,            dyt,              &
                                  tx(:,:,nt),     ty(:,:,nt))    

            do ij = 1, icells     ! ice is present
               i = indxi(ij)
               j = indxj(ij)
               tc(i,j,nt) = tm(i,j,nt)                    &
                          - tx(i,j,nt) * mtxav(i,j,nt1)   &
                          - ty(i,j,nt) * mtyav(i,j,nt1)
            enddo               ! ij

         elseif (tracer_type(nt)==3) then  ! upwind approx; gradient = 0

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               tc(i,j,nt) = tm(i,j,nt)
!               tx(i,j,nt) = c0   ! already initialized to 0.
!               ty(i,j,nt) = c0
            enddo               ! ij

         endif                  ! tracer_type
       enddo                    ! ntrace

      endif                     ! present (tm)

      end subroutine construct_fields

!=======================================================================
!
!BOP
!
! !IROUTINE: limited_gradient - limited gradient of a scalar field
!
! !INTERFACE:
!
      subroutine limited_gradient (nx_block, ny_block,   &
                                   ilo, ihi, jlo, jhi,   &
                                   nghost,               &
                                   phi,      phimask,    &
                                   cnx,      cny,        &
                                   HTN,      HTE,        &
                                   dxt,      dyt,        &
                                   gx,       gy)
!
! !DESCRIPTION:
!
! Compute a limited gradient of the scalar field phi in scaled coordinates.
! "Limited" means that we do not create new extrema in phi.  For
! instance, field values at the cell corners can neither exceed the
! maximum of phi(i,j) in the cell and its eight neighbors, nor fall
! below the minimum.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
          nx_block, ny_block,&! block dimensions
          ilo,ihi,jlo,jhi ,&! beginning and end of physical domain
          nghost              ! number of ghost cells

      real (kind=dbl_kind), dimension (nx_block,ny_block),   &
           intent (in) ::   &
          phi    ,&! input tracer field (mean values in each grid cell)
          cnx    ,&! x-coordinate of phi relative to geometric center of cell
          cny    ,&! y-coordinate of phi relative to geometric center of cell
          dxt    ,&! grid cell width (m)
          dyt    ,&! grid cell height (m)
          phimask ,&
          ! phimask(i,j) = 1 if phi(i,j) has physical meaning, = 0 otherwise.
          ! For instance, aice has no physical meaning in land cells,
          ! and hice no physical meaning where aice = 0.
          HTN    ,&! length of northern edge of T-cell (m)
          HTE      ! length of eastern edge of T-cell (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block),   &
          intent(out) ::   &
          gx     ,&! limited x-direction gradient
          gy       ! limited y-direction gradient
!
!EOP
!
      integer (kind=int_kind) ::   &
          i, j, ij        ,&! standard indices
          icells            ! number of cells to limit

      integer (kind=int_kind), dimension(nx_block*ny_block) ::   &
          indxi, indxj   ! combined i/j horizontal indices

      real (kind=dbl_kind) ::   &
          phi_nw, phi_n, phi_ne ,&! values of phi in 8 neighbor cells
          phi_w,         phi_e  ,&
          phi_sw, phi_s, phi_se ,&
          qmn, qmx     ,&! min and max value of phi within grid cell
          pmn, pmx     ,&! min and max value of phi among neighbor cells
          w1, w2, w3, w4 ! work variables

      real (kind=dbl_kind) ::   &
          gxtmp, gytmp   ! temporary term for x- and y- limited gradient

      gx(:,:) = c0
      gy(:,:) = c0

      ! For nghost = 1, loop over physical cells and update ghost cells later
      ! For nghost = 2, loop over a layer of ghost cells and skip the update

      icells = 0
      do j = jlo-nghost+1, jhi+nghost-1
      do i = ilo-nghost+1, ihi+nghost-1
         if (phimask(i,j) > puny) then

!jw            icells = icells + 1
!jw            indxi(icells) = i
!jw            indxj(icells) = j
!jw         endif                  ! phimask > puny
!jw      enddo
!jw      enddo

!jw      do ij = 1, icells
!jw         i = indxi(ij)
!jw         j = indxj(ij)

         ! Store values of phi in the 8 neighbor cells.
         ! Note: phimask = 1. or 0.  If phimask = 1., use the true value;
         !  if phimask = 0., use the home cell value so that non-physical
         !  values of phi do not contribute to the gradient.
         phi_nw = phimask(i-1,j+1) * phi(i-1,j+1)   &
            + (c1-phimask(i-1,j+1))* phi(i,j)
         phi_n  = phimask(i,j+1)   * phi(i,j+1)   &
            + (c1-phimask(i,j+1))  * phi(i,j)
         phi_ne = phimask(i+1,j+1) * phi(i+1,j+1)   &
            + (c1-phimask(i+1,j+1))* phi(i,j)
         phi_w  = phimask(i-1,j)   * phi(i-1,j)   &
            + (c1-phimask(i-1,j))  * phi(i,j)
         phi_e  = phimask(i+1,j)   * phi(i+1,j)   &
            + (c1-phimask(i+1,j))  * phi(i,j)
         phi_sw = phimask(i-1,j-1) * phi(i-1,j-1)   &
            + (c1-phimask(i-1,j-1))* phi(i,j)
         phi_s  = phimask(i,j-1)   * phi(i,j-1)   &
            + (c1-phimask(i,j-1))  * phi(i,j)
         phi_se = phimask(i+1,j-1) * phi(i+1,j-1)   &
            + (c1-phimask(i+1,j-1))* phi(i,j)

         ! unlimited gradient components
         ! (factors of two cancel out)

         gxtmp = (phi_e - phi(i,j)) / (dxt(i,j)   + dxt(i+1,j))   &
               + (phi(i,j) - phi_w) / (dxt(i-1,j) + dxt(i,j)  )
         gytmp = (phi_n - phi(i,j)) / (dyt(i,j)   + dyt(i,j+1))   &
               + (phi(i,j) - phi_s) / (dyt(i,j-1) + dyt(i,j)  )

         ! minimum and maximum among the nine local cells
         pmn = min (phi_nw, phi_n,  phi_ne, phi_w, phi(i,j),   &
                    phi_e,  phi_sw, phi_s,  phi_se)
         pmx = max (phi_nw, phi_n,  phi_ne, phi_w, phi(i,j),   &
                    phi_e,  phi_sw, phi_s,  phi_se)

         pmn = pmn - phi(i,j)
         pmx = pmx - phi(i,j)

         ! minimum and maximum deviation of phi within the cell

         w1  =  (p5*HTN(i,j)   - cnx(i,j)) * gxtmp   &
              + (p5*HTE(i,j)   - cny(i,j)) * gytmp
         w2  =  (p5*HTN(i,j-1) - cnx(i,j)) * gxtmp   &
              - (p5*HTE(i,j)   + cny(i,j)) * gytmp
         w3  = -(p5*HTN(i,j-1) + cnx(i,j)) * gxtmp   &
              - (p5*HTE(i-1,j) + cny(i,j)) * gytmp
         w4  =  (p5*HTE(i-1,j) - cny(i,j)) * gytmp   &
              - (p5*HTN(i,j)   + cnx(i,j)) * gxtmp

         qmn = min (w1, w2, w3, w4)
         qmx = max (w1, w2, w3, w4)

         ! Watch for underflows here

         ! the limiting coefficient
         if (abs(qmn) > 1.0e-300_dbl_kind) then ! 'abs(qmn) > puny' not sufficient
            w1 = max(c0, pmn/qmn)
         else
            w1 = c1
         endif

         if (abs(qmx) > 1.0e-300_dbl_kind) then
            w2 = max(c0, pmx/qmx)
         else
            w2 = c1
         endif

         w1 = min(c1, w1, w2)

         ! Limit the gradient components
         gx(i,j) = w1 * gxtmp
         gy(i,j) = w1 * gytmp

!jw      enddo                     ! ij
          endif
       enddo
      enddo

      end subroutine limited_gradient

!=======================================================================
!BOP
!
! !IROUTINE: departure_points - compute departure points of trajectories
!
! !INTERFACE:
!
      subroutine departure_points (nx_block,   ny_block,   &
                                   ilo, ihi,   jlo, jhi,   &
                                   nghost,     dt,   &
                                   uvel,       vvel,    &
                                   dxu,        dyu,     &
                                   HTN,        HTE,     &
                                   dpx,        dpy,     &
                                   l_dp_midpt, l_stop,   &
                                   istop,      jstop)
!
! !DESCRIPTION:
!
! Given velocity fields on cell corners, compute departure points
! of back trajectories in nondimensional coordinates.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         nx_block, ny_block,&! block dimensions
         ilo,ihi,jlo,jhi,   &! beginning and end of physical domain
         nghost              ! number of ghost cells

      real (kind=dbl_kind), intent(in) ::   &
         dt               ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::   &
         uvel           ,&! x-component of velocity (m/s)
         vvel           ,&! y-component of velocity (m/s)
         dxu            ,&! E-W dimensions of U-cell (m)
         dyu            ,&! N-S dimensions of U-cell (m)
         HTN            ,&! length of north face of T-cell (m) 
         HTE              ! length of east face of T-cell (m) 

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) ::   &
         dpx           ,&! coordinates of departure points (m)
         dpy             ! coordinates of departure points (m)

      logical (kind=log_kind), intent(in) ::   &
         l_dp_midpt          ! if true, find departure points using
                             ! corrected midpoint velocity

      logical (kind=log_kind), intent(inout) ::   &
         l_stop       ! if true, abort on return

      integer (kind=int_kind), intent(inout) ::   &
         istop, jstop     ! indices of grid cell where model aborts 
!
!EOP
!
      integer (kind=int_kind) ::   &
         i, j, i2, j2     ! horizontal indices

      real (kind=dbl_kind) ::                  &
         mpx,  mpy      ,&! coordinates of midpoint of back trajectory,
                          ! relative to cell corner
         mpxt, mpyt     ,&! midpoint coordinates relative to cell center
         ump,  vmp        ! corrected velocity at midpoint

    !-------------------------------------------------------------------
    ! Estimate departure points.
    ! This estimate is 1st-order accurate in time; improve accuracy by
    !  using midpoint approximation (to add later).
    ! For nghost = 1, loop over physical cells and update ghost cells later.
    ! For nghost = 2, loop over a layer of ghost cells and skip update.
    !-------------------------------------------------------------------

      dpx(:,:) = c0
      dpy(:,:) = c0

      do j = jlo-nghost+1, jhi+nghost-1
      do i = ilo-nghost+1, ihi+nghost-1

         dpx(i,j) = -dt*uvel(i,j)
         dpy(i,j) = -dt*vvel(i,j)

         ! Check for values out of bounds (more than one grid cell away)
         if (dpx(i,j) < -HTN(i,j) .or. dpx(i,j) > HTN(i+1,j) .or.   &
             dpy(i,j) < -HTE(i,j) .or. dpy(i,j) > HTE(i,j+1)) then
            l_stop = .true.
            istop = i
            jstop = j
         endif

      enddo
      enddo

      if (l_stop) then
         i = istop
         j = jstop
         write (nu_diag,*) ' '
         write (nu_diag,*)   &
                    'Warning: Departure points out of bounds in remap'
         write (nu_diag,*) 'my_task, i, j =', my_task, i, j
         write (nu_diag,*) 'dpx, dpy =', dpx(i,j), dpy(i,j)
         write (nu_diag,*) 'HTN(i,j), HTN(i+1,j) =', HTN(i,j), HTN(i+1,j)
         write (nu_diag,*) 'HTE(i,j), HTE(i,j+1) =', HTE(i,j), HTE(i,j+1)
         return
      endif

      if (l_dp_midpt) then ! find dep pts using corrected midpt velocity 

       do j = jlo-nghost+1, jhi+nghost-1
       do i = ilo-nghost+1, ihi+nghost-1
         if (uvel(i,j)/=c0 .or. vvel(i,j)/=c0) then
 
    !-------------------------------------------------------------------
    ! Scale departure points to coordinate system in which grid cells
    ! have sides of unit length.
    !-------------------------------------------------------------------

            dpx(i,j) = dpx(i,j) / dxu(i,j)
            dpy(i,j) = dpy(i,j) / dyu(i,j)

    !-------------------------------------------------------------------
    ! Estimate midpoint of backward trajectory relative to corner (i,j).
    !-------------------------------------------------------------------

            mpx = p5 * dpx(i,j)
            mpy = p5 * dpy(i,j)
 
    !-------------------------------------------------------------------
    ! Determine the indices (i2,j2) of the cell where the trajectory lies.
    ! Compute the coordinates of the midpoint of the backward trajectory
    !  relative to the cell center in a stretch coordinate system
    !  with vertices at (1/2, 1/2), (1/2, -1/2), etc.
    !-------------------------------------------------------------------

            if (mpx >= c0 .and. mpy >= c0) then    ! cell (i+1,j+1)
               i2 = i+1
               j2 = j+1
               mpxt = mpx - p5
               mpyt = mpy - p5
            elseif (mpx < c0 .and. mpy < c0) then  ! cell (i,j)
               i2 = i
               j2 = j
               mpxt = mpx + p5
               mpyt = mpy + p5
            elseif (mpx >= c0 .and. mpy < c0) then ! cell (i+1,j)
               i2 = i+1
               j2 = j
               mpxt = mpx - p5
               mpyt = mpy + p5
            elseif (mpx < c0 .and. mpy >= c0) then ! cell (i,j+1)
               i2 = i
               j2 = j+1
               mpxt = mpx + p5
               mpyt = mpy - p5
            endif
            
    !-------------------------------------------------------------------
    ! Using a bilinear approximation, estimate the velocity at the
    ! trajectory midpoint in the (i2,j2) reference frame.
    !-------------------------------------------------------------------
 
            ump = uvel(i2-1,j2-1)*(mpxt-p5)*(mpyt-p5)     &
                - uvel(i2,  j2-1)*(mpxt+p5)*(mpyt-p5)     &
                + uvel(i2,  j2  )*(mpxt+p5)*(mpyt+p5)     &  
                - uvel(i2-1,j2  )*(mpxt-p5)*(mpyt+p5)
 
            vmp = vvel(i2-1,j2-1)*(mpxt-p5)*(mpyt-p5)     &
                - vvel(i2,  j2-1)*(mpxt+p5)*(mpyt-p5)     &
                + vvel(i2,  j2  )*(mpxt+p5)*(mpyt+p5)     &
                - vvel(i2-1,j2  )*(mpxt-p5)*(mpyt+p5)
 
    !-------------------------------------------------------------------
    ! Use the midpoint velocity to estimate the coordinates of the
    !  departure point relative to corner (i,j).
    !-------------------------------------------------------------------
 
            dpx(i,j) = -dt * ump
            dpy(i,j) = -dt * vmp
 
         endif               ! nonzero velocity

       enddo                 ! i
       enddo                 ! j
 
      endif                  ! l_dp_midpt

      end subroutine departure_points

!=======================================================================
!
!BOP
!
! !IROUTINE: locate_triangles - triangle info for cell edges
!
! !INTERFACE:
!
      subroutine locate_triangles (nx_block,     ny_block,   &
                                   ilo, ihi,     jlo, jhi,   &
                                   nghost,       edge,       &
                                   icells,                   &
                                   indxi,        indxj,      &
                                   dpx,          dpy,        &
                                   dxu,          dyu,        &
                                   xp,           yp,         &
                                   iflux,        jflux,      &
                                   triarea,                  &
                                   l_fixed_area, edgearea)
!

! !DESCRIPTION:
!
! Compute areas and vertices of transport triangles for north or
!  east cell edges.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         nx_block, ny_block,&! block dimensions
         ilo,ihi,jlo,jhi   ,&! beginning and end of physical domain
         nghost              ! number of ghost cells

      character (len=char_len), intent(in) ::   &
         edge             ! 'north' or 'east'

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::  &
         dpx            ,&! x coordinates of departure points at cell corners
         dpy            ,&! y coordinates of departure points at cell corners
         dxu            ,&! E-W dimension of U-cell (m)
         dyu              ! N-S dimension of U-cell (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:nvert,ngroups),   &
         intent(out) ::   &
         xp, yp           ! coordinates of triangle vertices

      real (kind=dbl_kind), dimension (nx_block,ny_block,ngroups),   &
           intent(out) ::   &
         triarea          ! area of departure triangle

      integer (kind=int_kind), dimension (nx_block,ny_block,ngroups),    &
         intent(out) ::   &
         iflux          ,&! i index of cell contributing transport
         jflux            ! j index of cell contributing transport

      integer (kind=int_kind), dimension (ngroups), intent(out) ::   &
         icells           ! number of cells where triarea > puny

      integer (kind=int_kind), dimension (nx_block*ny_block,ngroups), &
         intent(out) ::                                               &
         indxi          ,&! compressed index in i-direction
         indxj            ! compressed index in j-direction

      logical, intent(in) ::   &
         l_fixed_area     ! if true, the area of each departure region is
                          !  passed in as edgearea
                          ! if false, edgearea if determined internally
                          !  and is passed out
                          
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) ::   &
         edgearea         ! area of departure region for each edge
                          ! edgearea > 0 for eastward/northward flow
!
!EOP
!
      integer (kind=int_kind) ::   &
         i, j, ij, ic   ,&! horizontal indices
         ib, ie, jb, je ,&! limits for loops over edges
         ng, nv         ,&! triangle indices
         ishift, jshift   ! differences between neighbor cells

      integer (kind=int_kind) ::   &
         icellsd          ! number of cells where departure area > 0.

      integer (kind=int_kind), dimension (nx_block*ny_block) ::  &
         indxid         ,&! compressed index in i-direction
         indxjd           ! compressed index in j-direction

      real (kind=dbl_kind), dimension(nx_block,ny_block) ::   &
         dx, dy         ,&! scaled departure points
         areafac_c      ,&! area scale factor at center of edge
         areafac_l      ,&! area scale factor at left corner
         areafac_r        ! area scale factor at right corner

      real (kind=dbl_kind) ::   &
         xcl, ycl       ,&! coordinates of left corner point
                          ! (relative to midpoint of edge)
         xdl, ydl       ,&! left departure point
         xil, yil       ,&! left intersection point
         xcr, ycr       ,&! right corner point
         xdr, ydr       ,&! right departure point
         xir, yir       ,&! right intersection point
         xic, yic       ,&! x-axis intersection point
         xicl, yicl     ,&! left-hand x-axis intersection point
         xicr, yicr     ,&! right-hand x-axis intersection point
         xdm, ydm       ,&! midpoint of segment connecting DL and DR;
                          ! shifted if l_fixed_area = T
         dxc            ,&! xcr - xcl
         dxd            ,&! xdr - xdl
         md             ,&! slope of line connecting DL and DR
         mdl            ,&! slope of line connecting DL and DM
         mdr            ,&! slope of line connecting DR and DM
         ishift_tl, jshift_tl ,&! i,j indices of TL cell relative to edge
         ishift_bl, jshift_bl ,&! i,j indices of BL cell relative to edge
         ishift_tr, jshift_tr ,&! i,j indices of TR cell relative to edge
         ishift_br, jshift_br ,&! i,j indices of BR cell relative to edge
         ishift_tc, jshift_tc ,&! i,j indices of TC cell relative to edge
         ishift_bc, jshift_bc ,&! i,j indices of BC cell relative to edge
         area1, area2         ,&! temporary triangle areas
         area3, area4         ,&! 
         area_c               ,&! center polygon area
         w1, w2                 ! work variables

      real (kind=dbl_kind), dimension (nx_block,ny_block,ngroups) ::   &
         areafact         ! = 1 for positive flux, -1 for negative

      real (kind=dbl_kind), dimension(nx_block,ny_block) ::   &
         areasum          ! sum of triangle areas for a given edge
      
    !-------------------------------------------------------------------
    ! Triangle notation:
    ! For each edge, there are 20 triangles that can contribute,
    ! but many of these are mutually exclusive.  It turns out that
    ! at most 5 triangles can contribute to transport integrals at once.
    !
    ! See Figure 3 in DB for pictures of these triangles.
    ! See Table 1 in DB for logical conditions.
    !
    ! For the north edge, DB refer to these triangles as:
    ! (1) NW, NW1, W, W2
    ! (2) NE, NE1, E, E2
    ! (3) NW2, W1, NE2, E1
    ! (4) H1a, H1b, N1a, N1b
    ! (5) H2a, H2b, N2a, N2b
    !
    ! For the east edge, DB refer to these triangles as:
    ! (1) NE, NE1, N, N2
    ! (2) SE, SE1, S, S2
    ! (3) NE2, N1, SE2, S1
    ! (4) H1a, H1b, E1a, E2b
    ! (5) H2a, H2b, E2a, E2b
    !
    ! The code below works for either north or east edges.
    ! The respective triangle labels are:
    ! (1) TL,  TL1, BL,  BL2
    ! (2) TR,  TR1, BR,  BR2
    ! (3) TL2, BL1, TR2, BR1
    ! (4) BC1a, BC1b, TC1a, TC2b
    ! (5) BC2a, BC2b, TC2a, TC2b
    ! 
    ! where the cell labels are:
    ! 
    !          |        |
    !     TL   |   TC   |   TR     (top left, center, right)
    !          |        |
    !   ------------------------
    !          |        |
    !     BL   |   BC   |   BR     (bottom left, center, right)
    !          |        |
    !
    ! and the transport is across the edge between cells TC and TB.
    !
    ! Departure points are scaled to a local coordinate system
    !  whose origin is at the midpoint of the edge.
    ! In this coordinate system, the lefthand corner CL = (-0.5,0)
    !  and the righthand corner CR = (0.5, 0).
    !-------------------------------------------------------------------
  
    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      areafac_c(:,:) = c0
      areafac_l(:,:) = c0
      areafac_r(:,:) = c0
      do ng = 1, ngroups
         do j = 1, ny_block
         do i = 1, nx_block
            triarea (i,j,ng) = c0
            areafact(i,j,ng) = c0
            iflux   (i,j,ng) = i
            jflux   (i,j,ng) = j
         enddo
         enddo
         do nv = 0, nvert
            do j = 1, ny_block
            do i = 1, nx_block
               xp(i,j,nv,ng) = c0
               yp(i,j,nv,ng) = c0
            enddo
            enddo
         enddo
      enddo

      if (trim(edge) == 'north') then

         ! loop size

         ib = ilo
         ie = ihi 
         jb = jlo - nghost            ! lowest j index is a ghost cell
         je = jhi

         ! index shifts for neighbor cells

         ishift_tl = -1
         jshift_tl =  1
         ishift_bl = -1
         jshift_bl =  0
         ishift_tr =  1
         jshift_tr =  1
         ishift_br =  1
         jshift_br =  0
         ishift_tc =  0
         jshift_tc =  1
         ishift_bc =  0
         jshift_bc =  0

         ! area scale factor

         do j = jb, je
         do i = ib, ie
            areafac_l(i,j) = dxu(i-1,j)*dyu(i-1,j) 
            areafac_r(i,j) = dxu(i,j)*dyu(i,j) 
            areafac_c(i,j) = p5*(areafac_l(i,j) + areafac_r(i,j))
         enddo
         enddo

      else                      ! east edge

         ! loop size

         ib = ilo - nghost            ! lowest i index is a ghost cell
         ie = ihi
         jb = jlo
         je = jhi

         ! index shifts for neighbor cells

         ishift_tl =  1
         jshift_tl =  1
         ishift_bl =  0
         jshift_bl =  1
         ishift_tr =  1
         jshift_tr = -1
         ishift_br =  0
         jshift_br = -1
         ishift_tc =  1
         jshift_tc =  0
         ishift_bc =  0
         jshift_bc =  0

         ! area scale factors

         do j = jb, je
         do i = ib, ie
            areafac_l(i,j) = dxu(i,j)*dyu(i,j) 
            areafac_r(i,j) = dxu(i,j-1)*dyu(i,j-1)
            areafac_c(i,j) = p5 * (areafac_l(i,j) + areafac_r(i,j))
         enddo
         enddo

      endif

    !-------------------------------------------------------------------
    ! Compute mask for edges with nonzero departure areas
    !-------------------------------------------------------------------

      if (l_fixed_area) then
         icellsd = 0
         do j = jb, je
         do i = ib, ie
            if (edgearea(i,j) /= c0) then
               icellsd = icellsd + 1
               indxid(icellsd) = i
               indxjd(icellsd) = j
            endif
         enddo
         enddo
      else
         icellsd = 0
         if (trim(edge) == 'north') then
            do j = jb, je
            do i = ib, ie
               if (dpx(i-1,j)/=c0 .or. dpy(i-1,j)/=c0   &
                                  .or.                  &
                     dpx(i,j)/=c0 .or.   dpy(i,j)/=c0) then
                  icellsd = icellsd + 1
                  indxid(icellsd) = i
                  indxjd(icellsd) = j
               endif
            enddo
            enddo
         else       ! east edge
            do j = jb, je
            do i = ib, ie
               if (dpx(i,j-1)/=c0 .or. dpy(i,j-1)/=c0   &
                                  .or.                  &
                     dpx(i,j)/=c0 .or.   dpy(i,j)/=c0) then
                  icellsd = icellsd + 1
                  indxid(icellsd) = i
                  indxjd(icellsd) = j
               endif
            enddo
            enddo
         endif       ! edge = north/east
      endif          ! l_fixed_area

    !-------------------------------------------------------------------
    ! Scale the departure points
    !-------------------------------------------------------------------

      do j = 1, je
      do i = 1, ie
         dx(i,j) = dpx(i,j) / dxu(i,j)
         dy(i,j) = dpy(i,j) / dyu(i,j)
      enddo
      enddo

    !-------------------------------------------------------------------
    ! Compute departure regions, divide into triangles, and locate
    !  vertices of each triangle.
    ! Work in a nondimensional coordinate system in which lengths are
    !  scaled by the local metric coefficients (dxu and dyu).
    ! Note: The do loop includes north faces of the j = 1 ghost cells
    !       when edge = 'north'.  The loop includes east faces of i = 1
    !       ghost cells when edge = 'east'.
    !-------------------------------------------------------------------

      do ij = 1, icellsd
         i = indxid(ij)
         j = indxjd(ij)
  
         xcl = -p5
         ycl =  c0

         xcr =  p5
         ycr =  c0

         ! Departure points

         if (trim(edge) == 'north') then ! north edge
            xdl = xcl + dx(i-1,j)
            ydl = ycl + dy(i-1,j)
            xdr = xcr + dx(i,j)
            ydr = ycr + dy(i,j)
         else                   ! east edge; rotate trajectory by pi/2
            xdl = xcl - dy(i,j)
            ydl = ycl + dx(i,j)
            xdr = xcr - dy(i,j-1)
            ydr = ycr + dx(i,j-1)
         endif

         xdm = p5 * (xdr + xdl)
         ydm = p5 * (ydr + ydl)

         ! Intersection points

         xil = xcl
         yil = (xcl*(ydm-ydl) + xdm*ydl - xdl*ydm) / (xdm - xdl)
         
         xir = xcr
         yir = (xcr*(ydr-ydm) - xdm*ydr + xdr*ydm) / (xdr - xdm) 
         
         md = (ydr - ydl) / (xdr - xdl)
         
         if (abs(md) > puny) then
            xic = xdl - ydl/md
         else
            xic = c0
         endif
         yic = c0

         xicl = xic
         yicl = yic
         xicr = xic
         yicr = yic

    !-------------------------------------------------------------------
    ! Locate triangles in TL cell (NW for north edge, NE for east edge)
    ! and BL cell (W for north edge, N for east edge).
    !-------------------------------------------------------------------

         if (yil > c0 .and. xdl < xcl .and. ydl >= c0) then

         ! TL (group 1)

            ng = 1
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xil
            yp    (i,j,2,ng) = yil
            xp    (i,j,3,ng) = xdl
            yp    (i,j,3,ng) = ydl
            iflux   (i,j,ng) = i + ishift_tl
            jflux   (i,j,ng) = j + jshift_tl
            areafact(i,j,ng) = -areafac_l(i,j)

         elseif (yil < c0 .and. xdl < xcl .and. ydl < c0) then

         ! BL (group 1)

            ng = 1
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xdl
            yp    (i,j,2,ng) = ydl
            xp    (i,j,3,ng) = xil
            yp    (i,j,3,ng) = yil
            iflux   (i,j,ng) = i + ishift_bl
            jflux   (i,j,ng) = j + jshift_bl
            areafact(i,j,ng) = areafac_l(i,j)

         elseif (yil < c0 .and. xdl < xcl .and. ydl >= c0) then

         ! TL1 (group 1)

            ng = 1
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xdl
            yp    (i,j,2,ng) = ydl
            xp    (i,j,3,ng) = xic
            yp    (i,j,3,ng) = yic
            iflux   (i,j,ng) = i + ishift_tl
            jflux   (i,j,ng) = j + jshift_tl
            areafact(i,j,ng) = areafac_l(i,j)

         ! BL1 (group 3)

            ng = 3
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xic
            yp    (i,j,2,ng) = yic
            xp    (i,j,3,ng) = xil
            yp    (i,j,3,ng) = yil
            iflux   (i,j,ng) = i + ishift_bl
            jflux   (i,j,ng) = j + jshift_bl
            areafact(i,j,ng) = areafac_l(i,j)

         elseif (yil > c0 .and. xdl < xcl .and. ydl < c0) then

         ! TL2 (group 3)

            ng = 3
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xil
            yp    (i,j,2,ng) = yil
            xp    (i,j,3,ng) = xic
            yp    (i,j,3,ng) = yic
            iflux   (i,j,ng) = i + ishift_tl
            jflux   (i,j,ng) = j + jshift_tl
            areafact(i,j,ng) = -areafac_l(i,j)

         ! BL2 (group 1)

            ng = 1
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xic
            yp    (i,j,2,ng) = yic
            xp    (i,j,3,ng) = xdl
            yp    (i,j,3,ng) = ydl
            iflux   (i,j,ng) = i + ishift_bl
            jflux   (i,j,ng) = j + jshift_bl
            areafact(i,j,ng) = -areafac_l(i,j)

         endif                  ! TL and BL triangles

    !-------------------------------------------------------------------
    ! Locate triangles in TR cell (NE for north edge, SE for east edge)
    ! and in BR cell (E for north edge, S for east edge).
    !-------------------------------------------------------------------

         if (yir > c0 .and. xdr >= xcr .and. ydr >= c0) then

         ! TR (group 2)

            ng = 2
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xdr
            yp    (i,j,2,ng) = ydr
            xp    (i,j,3,ng) = xir
            yp    (i,j,3,ng) = yir
            iflux   (i,j,ng) = i + ishift_tr
            jflux   (i,j,ng) = j + jshift_tr
            areafact(i,j,ng) = -areafac_r(i,j)

         elseif (yir < c0 .and. xdr >= xcr .and. ydr < c0) then

         ! BR (group 2)

            ng = 2
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xir
            yp    (i,j,2,ng) = yir
            xp    (i,j,3,ng) = xdr
            yp    (i,j,3,ng) = ydr
            iflux   (i,j,ng) = i + ishift_br
            jflux   (i,j,ng) = j + jshift_br
            areafact(i,j,ng) = areafac_r(i,j)

         elseif (yir < c0 .and. xdr >= xcr  .and. ydr >= c0) then 

         ! TR1 (group 2)

            ng = 2
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xic
            yp    (i,j,2,ng) = yic
            xp    (i,j,3,ng) = xdr
            yp    (i,j,3,ng) = ydr
            iflux   (i,j,ng) = i + ishift_tr
            jflux   (i,j,ng) = j + jshift_tr
            areafact(i,j,ng) = areafac_r(i,j)

         ! BR1 (group 3)

            ng = 3
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xir
            yp    (i,j,2,ng) = yir
            xp    (i,j,3,ng) = xic
            yp    (i,j,3,ng) = yic
            iflux   (i,j,ng) = i + ishift_br
            jflux   (i,j,ng) = j + jshift_br
            areafact(i,j,ng) = areafac_r(i,j)

         elseif (yir > c0 .and. xdr >= xcr .and. ydr < c0) then 

         ! TR2 (group 3)

            ng = 3
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xic
            yp    (i,j,2,ng) = yic
            xp    (i,j,3,ng) = xir
            yp    (i,j,3,ng) = yir
            iflux   (i,j,ng) = i + ishift_tr
            jflux   (i,j,ng) = j + jshift_tr
            areafact(i,j,ng) = -areafac_r(i,j)

         ! BR2 (group 2)

            ng = 2                     
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xdr
            yp    (i,j,2,ng) = ydr
            xp    (i,j,3,ng) = xic
            yp    (i,j,3,ng) = yic
            iflux   (i,j,ng) = i + ishift_br
            jflux   (i,j,ng) = j + jshift_br
            areafact(i,j,ng) = -areafac_r(i,j)

         endif                  ! TR and BR triangles

    !-------------------------------------------------------------------
    ! Redefine departure points if not located in central cells (TC or BC)
    !-------------------------------------------------------------------

         if (xdl < xcl) then
            xdl = xil
            ydl = yil
         endif

         if (xdr > xcr) then
            xdr = xir
            ydr = yir
         endif

    !-------------------------------------------------------------------
    ! For l_fixed_area = T, shift the midpoint so that the departure
    ! region has the prescribed area
    !-------------------------------------------------------------------

         if (l_fixed_area) then

            ! Sum the areas of the left and right triangles.
            ! Note that yp(i,j,1,ng) = 0 for all triangles, so we can
            !  drop those terms from the area formula.

            ng = 1
            area1 = p5 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
                            yp(i,j,3,ng)                   &
                         -  yp(i,j,2,ng) *                 &
                           (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
                         * areafact(i,j,ng) 

            ng = 2
            area2 = p5 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
                            yp(i,j,3,ng)                   &
                         -  yp(i,j,2,ng) *                 &
                           (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
                         * areafact(i,j,ng) 

            ng = 3
            area3 = p5 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
                            yp(i,j,3,ng)                   &
                         -  yp(i,j,2,ng) *                 &
                           (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
                         * areafact(i,j,ng) 

            !-----------------------------------------------------------
            ! Check whether the central triangles lie in one grid cell or two.
            ! If all are in one grid cell, then adjust the area of the central
            !  region so that the sum of all triangle areas is equal to the
            !  prescribed value.
            ! If two triangles are in one grid cell and one is in the other,
            !  then compute the area of the lone triangle using an area factor
            !  corresponding to the adjacent corner.  This is necessary to prevent
            !  negative masses in some rare cases on curved grids.  Then adjust
            !  the area of the remaining two-triangle region so that the sum of
            !  all triangle areas has the prescribed value.
            !-----------------------------------------------------------

            if (ydl*ydr >= c0) then   ! Both DPs lie on same side of x-axis

               ! compute required area of central departure region
               area_c  = edgearea(i,j) - area1 - area2 - area3

               ! shift midpoint so that the area of remaining triangles = area_c
               w1 = c2*area_c/areafac_c(i,j)    &
                    + (xdr-xcl)*ydl + (xcr-xdl)*ydr
               w2 = (xdr-xdl)**2 + (ydr-ydl)**2
               w1 = w1/w2
               xdm = xdm + (ydr - ydl) * w1
               ydm = ydm - (xdr - xdl) * w1

               ! compute left and right intersection points
               mdl = (ydm - ydl) / (xdm - xdl)
               mdr = (ydr - ydm) / (xdr - xdm)

               if (abs(mdl) > puny) then
                  xicl = xdl - ydl/mdl
               else
                  xicl = c0
               endif
               yicl = c0

               if (abs(mdr) > puny) then
                  xicr = xdr - ydr/mdr
               else
                  xicr = c0
               endif
               yicr = c0

            elseif (xic < c0) then  ! fix ICL = IC

               xicl = xic
               yicl = yic

               ! compute midpoint between ICL and DR
               xdm = p5 * (xdr + xicl)
               ydm = p5 *  ydr

               ! compute area of triangle adjacent to left corner 
               area4 = p5 * (xcl - xic) * ydl * areafac_l(i,j)
               area_c  = edgearea(i,j) - area1 - area2 - area3 - area4

               ! shift midpoint so that area of remaining triangles = area_c
               w1 = c2*area_c/areafac_c(i,j) + (xcr-xic)*ydr
               w2 = (xdr-xic)**2 + ydr**2
               w1 = w1/w2
               xdm = xdm + ydr*w1
               ydm = ydm - (xdr - xic) * w1

               ! compute ICR
               mdr = (ydr - ydm) / (xdr - xdm)
               if (abs(mdr) > puny) then
                  xicr = xdr - ydr/mdr
               else
                  xicr = c0
               endif
               yicr = c0

            elseif (xic >= c0) then  ! fix ICR = IR

               xicr = xic
               yicr = yic

               ! compute midpoint between ICR and DL 
               xdm = p5 * (xicr + xdl)
               ydm = p5 *  ydl

               area4 = p5 * (xic - xcr) * ydr * areafac_r(i,j)
               area_c  = edgearea(i,j) - area1 - area2 - area3 - area4

               ! shift midpoint so that area of remaining triangles = area_c
               w1 = c2*area_c/areafac_c(i,j) + (xic-xcl)*ydl
               w2 = (xic-xdl)**2 + ydl**2
               w1 = w1/w2
               xdm = xdm - ydl*w1
               ydm = ydm - (xic - xdl) * w1

               ! compute ICL

               mdl = (ydm - ydl) / (xdm - xdl)
               if (abs(mdl) > puny) then
                  xicl = xdl - ydl/mdl
               else
                  xicl = c0
               endif
               yicl = c0

            endif   ! ydl*ydr >= c0

         endif  ! l_fixed_area

    !-------------------------------------------------------------------
    ! Locate triangles in BC cell (H for both north and east edges) 
    ! and TC cell (N for north edge and E for east edge).
    !-------------------------------------------------------------------

    ! Start with cases where both DPs lie in the same grid cell

         if (ydl >= c0 .and. ydr >= c0 .and. ydm >= c0) then

         ! TC1a (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xcr
            yp    (i,j,2,ng) = ycr
            xp    (i,j,3,ng) = xdl
            yp    (i,j,3,ng) = ydl
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         ! TC2a (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xdr
            yp    (i,j,2,ng) = ydr
            xp    (i,j,3,ng) = xdl
            yp    (i,j,3,ng) = ydl
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         ! TC3a (group 6)
            ng = 6
            xp    (i,j,1,ng) = xdl
            yp    (i,j,1,ng) = ydl
            xp    (i,j,2,ng) = xdr
            yp    (i,j,2,ng) = ydr
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         elseif (ydl >= c0 .and. ydr >= c0 .and. ydm < c0) then  ! rare

         ! TC1b (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xicl
            yp    (i,j,2,ng) = yicl
            xp    (i,j,3,ng) = xdl
            yp    (i,j,3,ng) = ydl
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         ! TC2b (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xdr
            yp    (i,j,2,ng) = ydr
            xp    (i,j,3,ng) = xicr
            yp    (i,j,3,ng) = yicr
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         ! BC3b (group 6)

            ng = 6
            xp    (i,j,1,ng) = xicr
            yp    (i,j,1,ng) = yicr
            xp    (i,j,2,ng) = xicl
            yp    (i,j,2,ng) = yicl
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         elseif (ydl < c0 .and. ydr < c0 .and. ydm < c0) then

         ! BC1a (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xdl
            yp    (i,j,2,ng) = ydl
            xp    (i,j,3,ng) = xcr
            yp    (i,j,3,ng) = ycr
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         ! BC2a (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xdl
            yp    (i,j,2,ng) = ydl
            xp    (i,j,3,ng) = xdr
            yp    (i,j,3,ng) = ydr
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         ! BC3a (group 6)

            ng = 6
            xp    (i,j,1,ng) = xdl
            yp    (i,j,1,ng) = ydl
            xp    (i,j,2,ng) = xdm
            yp    (i,j,2,ng) = ydm
            xp    (i,j,3,ng) = xdr
            yp    (i,j,3,ng) = ydr
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         elseif (ydl < c0 .and. ydr < c0 .and. ydm >= c0) then  ! rare

         ! BC1b (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xdl
            yp    (i,j,2,ng) = ydl
            xp    (i,j,3,ng) = xicl
            yp    (i,j,3,ng) = yicl
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         ! BC2b (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xicr
            yp    (i,j,2,ng) = yicr
            xp    (i,j,3,ng) = xdr
            yp    (i,j,3,ng) = ydr
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         ! TC3b (group 6)

            ng = 6
            xp    (i,j,1,ng) = xicl
            yp    (i,j,1,ng) = yicl
            xp    (i,j,2,ng) = xicr
            yp    (i,j,2,ng) = yicr
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

    ! Now consider cases where the two DPs lie in different grid cells
    ! For these cases, one triangle is given the area factor associated
    !  with the adjacent corner, to avoid rare negative masses on curved grids.

         elseif (ydl >= c0 .and. ydr < c0 .and. xic >= c0  &
                                          .and. ydm >= c0) then

         ! TC1b (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xicr
            yp    (i,j,2,ng) = yicr
            xp    (i,j,3,ng) = xdl
            yp    (i,j,3,ng) = ydl
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         ! BC2b (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xicr
            yp    (i,j,2,ng) = yicr
            xp    (i,j,3,ng) = xdr
            yp    (i,j,3,ng) = ydr
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_r(i,j)

         ! TC3b (group 6)

            ng = 6
            xp    (i,j,1,ng) = xdl
            yp    (i,j,1,ng) = ydl
            xp    (i,j,2,ng) = xicr
            yp    (i,j,2,ng) = yicr
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         elseif (ydl >= c0 .and. ydr < c0 .and. xic >= c0  &
                                          .and. ydm < c0 ) then  ! less common

         ! TC1b (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xicl
            yp    (i,j,2,ng) = yicl
            xp    (i,j,3,ng) = xdl
            yp    (i,j,3,ng) = ydl
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         ! BC2b (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xicr
            yp    (i,j,2,ng) = yicr
            xp    (i,j,3,ng) = xdr
            yp    (i,j,3,ng) = ydr
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_r(i,j)

         ! BC3b (group 6)

            ng = 6
            xp    (i,j,1,ng) = xicr
            yp    (i,j,1,ng) = yicr
            xp    (i,j,2,ng) = xicl
            yp    (i,j,2,ng) = yicl
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         elseif (ydl >= c0 .and. ydr < c0 .and. xic < c0   &
                                          .and. ydm < c0) then

         ! TC1b (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xicl
            yp    (i,j,2,ng) = yicl
            xp    (i,j,3,ng) = xdl
            yp    (i,j,3,ng) = ydl
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_l(i,j)

         ! BC2b (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xicl
            yp    (i,j,2,ng) = yicl
            xp    (i,j,3,ng) = xdr
            yp    (i,j,3,ng) = ydr
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         ! BC3b (group 6)

            ng = 6
            xp    (i,j,1,ng) = xdr
            yp    (i,j,1,ng) = ydr
            xp    (i,j,2,ng) = xicl
            yp    (i,j,2,ng) = yicl
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         elseif (ydl >= c0 .and. ydr < c0 .and. xic <  c0  &
                                          .and. ydm >= c0) then  ! less common

         ! TC1b (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xicl
            yp    (i,j,2,ng) = yicl
            xp    (i,j,3,ng) = xdl
            yp    (i,j,3,ng) = ydl
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_l(i,j)

         ! BC2b (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xicr
            yp    (i,j,2,ng) = yicr
            xp    (i,j,3,ng) = xdr
            yp    (i,j,3,ng) = ydr
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         ! TC3b (group 6)

            ng = 6
            xp    (i,j,1,ng) = xicl
            yp    (i,j,1,ng) = yicl
            xp    (i,j,2,ng) = xicr
            yp    (i,j,2,ng) = yicr
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         elseif (ydl < c0 .and. ydr >= c0 .and. xic <  c0  &
                                          .and. ydm >= c0) then

         ! BC1b (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xdl
            yp    (i,j,2,ng) = ydl
            xp    (i,j,3,ng) = xicl
            yp    (i,j,3,ng) = yicl
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_l(i,j)

         ! TC2b (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xdr
            yp    (i,j,2,ng) = ydr
            xp    (i,j,3,ng) = xicl
            yp    (i,j,3,ng) = yicl
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         ! TC3b (group 6)

            ng = 6
            xp    (i,j,1,ng) = xicl
            yp    (i,j,1,ng) = yicl
            xp    (i,j,2,ng) = xdr
            yp    (i,j,2,ng) = ydr
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         elseif (ydl < c0 .and. ydr >= c0 .and. xic < c0  &
                                          .and. ydm < c0) then ! less common

         ! BC1b (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xdl
            yp    (i,j,2,ng) = ydl
            xp    (i,j,3,ng) = xicl
            yp    (i,j,3,ng) = yicl
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_l(i,j)

         ! TC2b (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xdr
            yp    (i,j,2,ng) = ydr
            xp    (i,j,3,ng) = xicr
            yp    (i,j,3,ng) = yicr
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         ! BC3b (group 6)

            ng = 6
            xp    (i,j,1,ng) = xicr
            yp    (i,j,1,ng) = yicr
            xp    (i,j,2,ng) = xicl
            yp    (i,j,2,ng) = yicl
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         elseif (ydl < c0 .and. ydr >= c0 .and. xic >= c0  &
                                          .and. ydm <  c0) then

         ! BC1b (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xdl
            yp    (i,j,2,ng) = ydl
            xp    (i,j,3,ng) = xicr
            yp    (i,j,3,ng) = yicr
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         ! TC2b (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xdr
            yp    (i,j,2,ng) = ydr
            xp    (i,j,3,ng) = xicr
            yp    (i,j,3,ng) = yicr
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_r(i,j)

         ! BC3b (group 6)

            ng = 6
            xp    (i,j,1,ng) = xicr
            yp    (i,j,1,ng) = yicr
            xp    (i,j,2,ng) = xdl
            yp    (i,j,2,ng) = ydl
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         elseif (ydl < c0 .and. ydr >= c0 .and. xic >= c0   &
                                          .and. ydm >= c0) then  ! less common

         ! BC1b (group 4)

            ng = 4
            xp    (i,j,1,ng) = xcl
            yp    (i,j,1,ng) = ycl
            xp    (i,j,2,ng) = xdl
            yp    (i,j,2,ng) = ydl
            xp    (i,j,3,ng) = xicl
            yp    (i,j,3,ng) = yicl
            iflux   (i,j,ng) = i + ishift_bc
            jflux   (i,j,ng) = j + jshift_bc
            areafact(i,j,ng) = areafac_c(i,j)

         ! TC2b (group 5)

            ng = 5
            xp    (i,j,1,ng) = xcr
            yp    (i,j,1,ng) = ycr
            xp    (i,j,2,ng) = xdr
            yp    (i,j,2,ng) = ydr
            xp    (i,j,3,ng) = xicr
            yp    (i,j,3,ng) = yicr
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_r(i,j)

         ! TC3b (group 6)

            ng = 6
            xp    (i,j,1,ng) = xicl
            yp    (i,j,1,ng) = yicl
            xp    (i,j,2,ng) = xicr
            yp    (i,j,2,ng) = yicr
            xp    (i,j,3,ng) = xdm
            yp    (i,j,3,ng) = ydm
            iflux   (i,j,ng) = i + ishift_tc
            jflux   (i,j,ng) = j + jshift_tc
            areafact(i,j,ng) = -areafac_c(i,j)

         endif                  ! TC and BC triangles

      enddo                     ! ij

    !-------------------------------------------------------------------
    ! Compute triangle areas with appropriate sign.
    ! These are found by computing the area in scaled coordinates and
    !  multiplying by a scale factor (areafact).
    ! Note that the scale factor is positive for fluxes out of the cell 
    !  and negative for fluxes into the cell.
    !
    ! Note: The triangle area formula below gives A >=0 iff the triangle
    !        points x1, x2, and x3 are taken in counterclockwise order.
    !       These points are defined above in such a way that the
    !        order is nearly always CCW.
    !       In rare cases, we may compute A < 0.  In this case,
    !        the quadrilateral departure area is equal to the 
    !        difference of two triangle areas instead of the sum.
    !        The fluxes work out correctly in the end.
    !
    ! Also compute the cumulative area transported across each edge.
    ! If l_fixed_area = T, this area is compared to edgearea as a bug check.
    ! If l_fixed_area = F, this area is passed as an output array.
    !-------------------------------------------------------------------

      areasum(:,:) = c0

      do ng = 1, ngroups
         icells(ng) = 0

         do ij = 1, icellsd
            i = indxid(ij)
            j = indxjd(ij)

            triarea(i,j,ng) = p5 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
                                     (yp(i,j,3,ng)-yp(i,j,1,ng))   &
                                   - (yp(i,j,2,ng)-yp(i,j,1,ng)) *   &
                                     (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
                                   * areafact(i,j,ng) 

            if (abs(triarea(i,j,ng)) < eps16*areafac_c(i,j)) then
               triarea(i,j,ng) = c0
            else
               icells(ng) = icells(ng) + 1 
               ic = icells(ng)
               indxi(ic,ng) = i
               indxj(ic,ng) = j
            endif

            areasum(i,j) = areasum(i,j) + triarea(i,j,ng)

         enddo                  ! ij
      enddo                     ! ng

      if (l_fixed_area) then
       if (bugcheck) then   ! set bugcheck = F to speed up code
         do ij = 1, icellsd
            i = indxid(ij)
            j = indxjd(ij)
            if (abs(areasum(i,j) - edgearea(i,j)) > eps13*areafac_c(i,j)) then
               print*, ''
               print*, 'Areas do not add up: m, i, j, edge =',   &
                        my_task, i, j, trim(edge)
               print*, 'edgearea =', edgearea(i,j)
               print*, 'areasum =', areasum(i,j)
               print*, 'areafac_c =', areafac_c(i,j)
               print*, ''
               print*, 'Triangle areas:'
               do ng = 1, ngroups   ! not vector friendly
                  if (abs(triarea(i,j,ng)) > eps16*abs(areafact(i,j,ng))) then
                     print*, ng, triarea(i,j,ng)
                  endif
               enddo
            endif
         enddo
       endif          ! bugcheck

      else            ! l_fixed_area = F
         do ij = 1, icellsd
            i = indxid(ij)
            j = indxjd(ij)
            edgearea(i,j) = areasum(i,j)
         enddo
      endif     ! l_fixed_area

    !-------------------------------------------------------------------
    ! Transform triangle vertices to a scaled coordinate system centered
    !  in the cell containing the triangle.
    !-------------------------------------------------------------------

      if (trim(edge) == 'north') then
         do ng = 1, ngroups
            do nv = 1, nvert
               do ij = 1, icells(ng)
                  i = indxi(ij,ng)
                  j = indxj(ij,ng)
                  ishift = iflux(i,j,ng) - i
                  jshift = jflux(i,j,ng) - j
                  xp(i,j,nv,ng) = xp(i,j,nv,ng) - c1*ishift
                  yp(i,j,nv,ng) = yp(i,j,nv,ng) + p5 - c1*jshift
               enddo            ! ij
            enddo               ! nv
         enddo                  ! ng
      else                      ! east edge
         do ng = 1, ngroups
            do nv = 1, nvert
               do ij = 1, icells(ng)
                  i = indxi(ij,ng)
                  j = indxj(ij,ng)
                  ishift = iflux(i,j,ng) - i
                  jshift = jflux(i,j,ng) - j
                  ! Note rotation of pi/2 here
                  w1 = xp(i,j,nv,ng)
                  xp(i,j,nv,ng) =  yp(i,j,nv,ng) + p5 - c1*ishift
                  yp(i,j,nv,ng) = -w1 - c1*jshift
               enddo            ! ij
            enddo               ! nv
         enddo                  ! ng
      endif

      if (bugcheck) then
         do ng = 1, ngroups
         do nv = 1, nvert
            do j = jb, je
            do i = ib, ie
               if (abs(triarea(i,j,ng)) > puny) then
                  if (abs(xp(i,j,nv,ng)) > p5+puny) then
                     print*, ''
                     print*, 'WARNING: xp =', xp(i,j,nv,ng)
                     print*, 'm, i, j, ng, nv =', my_task, i, j, ng, nv
!                     print*, 'yil,xdl,xcl,ydl=',yil,xdl,xcl,ydl
!                     print*, 'yir,xdr,xcr,ydr=',yir,xdr,xcr,ydr
!                     print*, 'ydm=',ydm
!                      stop
                  endif
                  if (abs(yp(i,j,nv,ng)) > p5+puny) then
                     print*, ''
                     print*, 'WARNING: yp =', yp(i,j,nv,ng)
                     print*, 'm, i, j, ng, nv =', my_task, i, j, ng, nv
                  endif
               endif   ! triarea
            enddo
            enddo
         enddo
         enddo
      endif  ! bugcheck

      end subroutine locate_triangles

!=======================================================================
!
!BOP
! !IROUTINE: triangle_coordinates - find coordinates of quadrature points
!
! !INTERFACE:
!
      subroutine triangle_coordinates (nx_block,       ny_block,  &
                                       integral_order, icells,    &
                                       indxi,          indxj,     &
                                       xp,             yp)
!
! !DESCRIPTION:
!
! For each triangle, find the coordinates of the quadrature points needed
!  to compute integrals of linear, quadratic, or cubic polynomials,
!  using formulas from A.H. Stroud, Approximate Calculation of Multiple
!  Integrals, Prentice-Hall, 1971.  (Section 8.8, formula 3.1.)
! Linear functions can be integrated exactly by evaluating the function 
!  at just one point (the midpoint).  Quadratic functions require
!  3 points, and cubics require 4 points.
! The default is cubic, but the code can be sped up slightly using 
!  linear or quadratic integrals, usually with little loss of accuracy.
!
! The formulas are as follows:
!
! I1 = integral of f(x,y)*dA
!    = A * f(x0,y0)
! where A is the traingle area and (x0,y0) is the midpoint.
!
! I2 = A * (f(x1,y1) + f(x2,y2) + f(x3,y3))
! where these three points are located halfway between the midpoint
! and the three vertics of the triangle.
!
! I3 = A * [ -9/16 *  f(x0,y0)
!           + 25/48 * (f(x1,y1) + f(x2,y2) + f(x3,y3))]
! where (x0,y0) is the midpoint, and the other three points are
! located 2/5 of the way from the midpoint to the three vertices.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
           nx_block, ny_block,&! block dimensions
           integral_order      ! polynomial order for quadrature integrals 

      integer (kind=int_kind), dimension (ngroups), intent(in) ::     &
           icells              ! number of cells where triarea > puny

      integer (kind=int_kind), dimension (nx_block*ny_block,ngroups),     &
           intent(in) ::     &
           indxi ,&! compressed index in i-direction
           indxj   ! compressed index in j-direction

      real (kind=dbl_kind), intent(inout),   &
           dimension (nx_block, ny_block, 0:nvert, ngroups) ::   &
           xp, yp          ! coordinates of triangle points
!
!EOP
!
      integer (kind=int_kind) ::   &
           i, j, ij          ,&! horizontal indices
           ng                  ! triangle index


      if (integral_order == 1) then ! linear (1-point formula)

         do ng = 1, ngroups
         do ij = 1, icells(ng)
            i = indxi(ij,ng)
            j = indxj(ij,ng)

            ! coordinates of midpoint
            xp(i,j,0,ng) = p333   &
                        * (xp(i,j,1,ng) + xp(i,j,2,ng) + xp(i,j,3,ng))
            yp(i,j,0,ng) = p333   &
                        * (yp(i,j,1,ng) + yp(i,j,2,ng) + yp(i,j,3,ng))

         enddo                  ! ij
         enddo                  ! ng

      elseif (integral_order == 2) then ! quadratic (3-point formula)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

         do ng = 1, ngroups
         do ij = 1, icells(ng)
            i = indxi(ij,ng)
            j = indxj(ij,ng)

            ! coordinates of midpoint
            xp(i,j,0,ng) = p333   &
                        * (xp(i,j,1,ng) + xp(i,j,2,ng) + xp(i,j,3,ng))
            yp(i,j,0,ng) = p333   &
                        * (yp(i,j,1,ng) + yp(i,j,2,ng) + yp(i,j,3,ng))

            ! coordinates of the 3 points needed for integrals

            xp(i,j,1,ng) = p5*xp(i,j,1,ng) + p5*xp(i,j,0,ng)
            yp(i,j,1,ng) = p5*yp(i,j,1,ng) + p5*yp(i,j,0,ng)

            xp(i,j,2,ng) = p5*xp(i,j,2,ng) + p5*xp(i,j,0,ng)
            yp(i,j,2,ng) = p5*yp(i,j,2,ng) + p5*yp(i,j,0,ng)

            xp(i,j,3,ng) = p5*xp(i,j,3,ng) + p5*xp(i,j,0,ng)
            yp(i,j,3,ng) = p5*yp(i,j,3,ng) + p5*yp(i,j,0,ng)

         enddo                  ! ij
         enddo                  ! ng

      else                      ! cubic (4-point formula)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ng = 1, ngroups
         do ij = 1, icells(ng)
            i = indxi(ij,ng)
            j = indxj(ij,ng)

            ! coordinates of midpoint
            xp(i,j,0,ng) = p333   &
                        * (xp(i,j,1,ng) + xp(i,j,2,ng) + xp(i,j,3,ng))
            yp(i,j,0,ng) = p333   &
                        * (yp(i,j,1,ng) + yp(i,j,2,ng) + yp(i,j,3,ng))

            ! coordinates of the other 3 points needed for integrals

            xp(i,j,1,ng) = p4*xp(i,j,1,ng) + p6*xp(i,j,0,ng)
            yp(i,j,1,ng) = p4*yp(i,j,1,ng) + p6*yp(i,j,0,ng)

            xp(i,j,2,ng) = p4*xp(i,j,2,ng) + p6*xp(i,j,0,ng)
            yp(i,j,2,ng) = p4*yp(i,j,2,ng) + p6*yp(i,j,0,ng)
               
            xp(i,j,3,ng) = p4*xp(i,j,3,ng) + p6*xp(i,j,0,ng)
            yp(i,j,3,ng) = p4*yp(i,j,3,ng) + p6*yp(i,j,0,ng)
               
         enddo                  ! ij
         enddo                  ! ng

      endif

      end subroutine triangle_coordinates

!=======================================================================
!
!BOP
!
! !IROUTINE: transport_integrals - compute transports across each edge
!
! !INTERFACE:
!
      subroutine transport_integrals (nx_block,       ny_block,    &
                                      ntrace,         icells,      &
                                      indxi,          indxj,       &
                                      tracer_type,    depend,      &
                                      integral_order, triarea,     &
                                      iflux,          jflux,       &
                                      xp,             yp,          &
                                      mc,             mx,          &
                                      my,             mflx,       &
                                      tc,             tx,          &
                                      ty,             mtflx)
!
! !DESCRIPTION:
!
! Compute the transports across each edge by integrating the mass
! and tracers over each departure triangle.
! Input variables have the same meanings as in the main subroutine.
! Repeated use of certain sums makes the calculation more efficient.
! Integral formulas are described in triangle_coordinates subroutine.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
           nx_block, ny_block  ,&! block dimensions
           ntrace              ,&! number of tracers in use
           integral_order   ! polynomial order for quadrature integrals 

      integer (kind=int_kind), dimension (ntrace), intent(in) ::     &
           tracer_type       ,&! = 1, 2, or 3 (see comments above)
           depend              ! tracer dependencies (see above)

      integer (kind=int_kind), dimension (ngroups), intent(in) ::     &
           icells           ! number of cells where triarea > puny

      integer (kind=int_kind), dimension (nx_block*ny_block,ngroups),     &
           intent(in) ::     &
           indxi ,&! compressed index in i-direction
           indxj   ! compressed index in j-direction

      real (kind=dbl_kind), intent(in),   &
           dimension (nx_block, ny_block, 0:nvert, ngroups) ::   &
           xp, yp           ! coordinates of triangle points

      real (kind=dbl_kind), intent(in),   &
           dimension (nx_block, ny_block, ngroups) ::   &
           triarea          ! triangle area

      integer (kind=int_kind), intent(in),   &
           dimension (nx_block, ny_block, ngroups) ::   &
           iflux     ,&
           jflux

      real (kind=dbl_kind), intent(in),   &
           dimension (nx_block, ny_block) ::   &
           mc, mx, my

      real (kind=dbl_kind), intent(out),   &
           dimension (nx_block, ny_block) ::   &
           mflx

      real (kind=dbl_kind), intent(in),   &
           dimension (nx_block, ny_block, max_ntrace), optional ::   &
           tc, tx, ty

      real (kind=dbl_kind), intent(out),   &
           dimension (nx_block, ny_block, max_ntrace), optional ::   &
           mtflx
!
!EOP
!
      integer (kind=int_kind) ::   &
           i, j, ij      ,&! horizontal indices of edge
           i2, j2        ,&! horizontal indices of cell contributing transport
           ng            ,&! triangle index
           nt, nt1       ,&! tracer indices
           ilo,ihi,jlo,jhi ! beginning and end of physical domain

      real (kind=dbl_kind) ::   &
           m0, m1, m2, m3         ,&! mass field at internal points
           w0, w1, w2, w3           ! work variables

      real (kind=dbl_kind), dimension (nx_block, ny_block) ::   &
           msum, mxsum, mysum     ,&! sum of mass, mass*x, and mass*y
           mxxsum, mxysum, myysum   ! sum of mass*x*x, mass*x*y, mass*y*y

      real (kind=dbl_kind), dimension (nx_block, ny_block, max_ntrace) ::   &
           mtsum            ,&! sum of mass*tracer
           mtxsum           ,&! sum of mass*tracer*x
           mtysum             ! sum of mass*tracer*y

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      mflx(:,:) = c0
      if (present(mtflx)) then
         do nt = 1, ntrace
            mtflx(:,:,nt) = c0
         enddo
      endif

    !-------------------------------------------------------------------
    ! Main loop
    !-------------------------------------------------------------------

      do ng = 1, ngroups

         if (integral_order == 1) then  ! linear (1-point formula)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells(ng)
               i = indxi(ij,ng)
               j = indxj(ij,ng)

               i2 = iflux(i,j,ng)
               j2 = jflux(i,j,ng)

               ! mass transports

               m0 = mc(i2,j2) + xp(i,j,0,ng)*mx(i2,j2)   &
                              + yp(i,j,0,ng)*my(i2,j2)
               msum(i,j) = m0

               mflx(i,j) = mflx(i,j) + triarea(i,j,ng)*msum(i,j)

               ! quantities needed for tracer transports
               mxsum(i,j)  =         m0*xp(i,j,0,ng) 
               mxxsum(i,j) = mxsum(i,j)*xp(i,j,0,ng) 
               mxysum(i,j) = mxsum(i,j)*yp(i,j,0,ng) 
               mysum(i,j)  =         m0*yp(i,j,0,ng) 
               myysum(i,j) = mysum(i,j)*yp(i,j,0,ng) 
            enddo               ! ij

         elseif (integral_order == 2) then  ! quadratic (3-point formula)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells(ng)
               i = indxi(ij,ng)
               j = indxj(ij,ng)

               i2 = iflux(i,j,ng)
               j2 = jflux(i,j,ng)

               ! mass transports
               ! Weighting factor of 1/3 is incorporated into the ice
               ! area terms m1, m2, and m3.
               m1 = p333 * (mc(i2,j2) + xp(i,j,1,ng)*mx(i2,j2)   &
                                      + yp(i,j,1,ng)*my(i2,j2))
               m2 = p333 * (mc(i2,j2) + xp(i,j,2,ng)*mx(i2,j2)   &
                                      + yp(i,j,2,ng)*my(i2,j2))
               m3 = p333 * (mc(i2,j2) + xp(i,j,3,ng)*mx(i2,j2)   &
                                      + yp(i,j,3,ng)*my(i2,j2))
               msum(i,j) = m1 + m2 + m3
               mflx(i,j) = mflx(i,j) + triarea(i,j,ng)*msum(i,j)

               ! quantities needed for mass_tracer transports
               w1 = m1 * xp(i,j,1,ng)
               w2 = m2 * xp(i,j,2,ng)
               w3 = m3 * xp(i,j,3,ng)

               mxsum(i,j) = w1 + w2 + w3

               mxxsum(i,j) = w1*xp(i,j,1,ng) + w2*xp(i,j,2,ng)   &
                           + w3*xp(i,j,3,ng) 

               mxysum(i,j) = w1*yp(i,j,1,ng) + w2*yp(i,j,2,ng)   &
                           + w3*yp(i,j,3,ng)

               w1 = m1 * yp(i,j,1,ng)
               w2 = m2 * yp(i,j,2,ng)
               w3 = m3 * yp(i,j,3,ng)

               mysum(i,j) = w1 + w2 + w3

               myysum(i,j) = w1*yp(i,j,1,ng) + w2*yp(i,j,2,ng)   &
                           + w3*yp(i,j,3,ng)
            enddo               ! ij

         else                   ! cubic (4-point formula)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells(ng)
               i = indxi(ij,ng)
               j = indxj(ij,ng)

               i2 = iflux(i,j,ng)
               j2 = jflux(i,j,ng)

               ! mass transports

               ! Weighting factors are incorporated into the
               ! terms m0, m1, m2, and m3.
               m0 = p5625m * (mc(i2,j2) + xp(i,j,0,ng)*mx(i2,j2)   &
                                        + yp(i,j,0,ng)*my(i2,j2))
               m1 = p52083 * (mc(i2,j2) + xp(i,j,1,ng)*mx(i2,j2)   &
                                        + yp(i,j,1,ng)*my(i2,j2))
               m2 = p52083 * (mc(i2,j2) + xp(i,j,2,ng)*mx(i2,j2)   &
                                        + yp(i,j,2,ng)*my(i2,j2))
               m3 = p52083 * (mc(i2,j2) + xp(i,j,3,ng)*mx(i2,j2)   &
                                        + yp(i,j,3,ng)*my(i2,j2))
               msum(i,j) = m0 + m1 + m2 + m3
               mflx(i,j) = mflx(i,j) + triarea(i,j,ng)*msum(i,j)

               ! quantities needed for tracer transports
               w0 = m0 * xp(i,j,0,ng)
               w1 = m1 * xp(i,j,1,ng)
               w2 = m2 * xp(i,j,2,ng)
               w3 = m3 * xp(i,j,3,ng)

               mxsum(i,j) = w0 + w1 + w2 + w3

               mxxsum(i,j) = w0*xp(i,j,0,ng) + w1*xp(i,j,1,ng)   &
                           + w2*xp(i,j,2,ng) + w3*xp(i,j,3,ng)

               mxysum(i,j) = w0*yp(i,j,0,ng) + w1*yp(i,j,1,ng)   &
                           + w2*yp(i,j,2,ng) + w3*yp(i,j,3,ng)

               w0 = m0 * yp(i,j,0,ng)
               w1 = m1 * yp(i,j,1,ng)
               w2 = m2 * yp(i,j,2,ng)
               w3 = m3 * yp(i,j,3,ng)

               mysum(i,j) = w0 + w1 + w2 + w3

               myysum(i,j) = w0*yp(i,j,0,ng) + w1*yp(i,j,1,ng)   &
                           + w2*yp(i,j,2,ng) + w3*yp(i,j,3,ng)

            enddo               ! ij

         endif                  ! integral_order

         ! mass * tracer transports

         if (present(mtflx)) then

            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, icells(ng)
                     i = indxi(ij,ng)
                     j = indxj(ij,ng)

                     i2 = iflux(i,j,ng)
                     j2 = jflux(i,j,ng)

                     mtsum(i,j,nt) =  msum(i,j) * tc(i2,j2,nt)   &
                                   + mxsum(i,j) * tx(i2,j2,nt)   &
                                   + mysum(i,j) * ty(i2,j2,nt)

                     mtflx(i,j,nt) = mtflx(i,j,nt)   &
                                 + triarea(i,j,ng) * mtsum(i,j,nt)

                     ! quantities needed for dependent tracers

                     mtxsum(i,j,nt) =  mxsum(i,j) * tc(i2,j2,nt)   &
                                    + mxxsum(i,j) * tx(i2,j2,nt)   &
                                    + mxysum(i,j) * ty(i2,j2,nt)

                     mtysum(i,j,nt) =  mysum(i,j) * tc(i2,j2,nt)   &
                                    + mxysum(i,j) * tx(i2,j2,nt)   &
                                    + myysum(i,j) * ty(i2,j2,nt)
                  enddo         ! ij

               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, icells(ng)
                     i = indxi(ij,ng)
                     j = indxj(ij,ng)

                     i2 = iflux(i,j,ng)
                     j2 = jflux(i,j,ng)

                     mtsum(i,j,nt) =  mtsum(i,j,nt1) * tc(i2,j2,nt)   &
                                   + mtxsum(i,j,nt1) * tx(i2,j2,nt)   &
                                   + mtysum(i,j,nt1) * ty(i2,j2,nt)

                     mtflx(i,j,nt) = mtflx(i,j,nt)   &
                                   + triarea(i,j,ng) * mtsum(i,j,nt)
                  enddo         ! ij


               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, icells(ng)
                     i = indxi(ij,ng)
                     j = indxj(ij,ng)

                     i2 = iflux(i,j,ng)
                     j2 = jflux(i,j,ng)

                     ! upwind approx (tx=ty=0) for type 3 tracers
                     mtsum(i,j,nt) =  mtsum(i,j,nt1) * tc(i2,j2,nt)

                     mtflx(i,j,nt) = mtflx(i,j,nt)   &
                                   + triarea(i,j,ng) * mtsum(i,j,nt)
                  enddo         ! ij

               endif            ! tracer type
            enddo               ! ntrace
         endif                  ! present(mtflx)
      enddo                     ! ng

      end subroutine transport_integrals

!=======================================================================
!
!BOP
!
! !IROUTINE: update_fields - compute new area and tracers
!
! !INTERFACE:
!
      subroutine update_fields (nx_block,    ny_block,   &
                                ilo, ihi,    jlo, jhi,   &
                                ntrace,                  &
                                tracer_type, depend,     &
                                tarear,      l_stop,     &
                                istop,       jstop,      &
                                mflxe,       mflxn,      &
                                mm,                      &
                                mtflxe,      mtflxn,     &
                                tm)
!
! !DESCRIPTION:
!
! Given transports through cell edges, compute new area and tracers.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::   &
         nx_block, ny_block,&! block dimensions
         ilo,ihi,jlo,jhi   ,&! beginning and end of physical domain
         ntrace              ! number of tracers in use

      integer (kind=int_kind), dimension (ntrace), intent(in) ::     &
         tracer_type       ,&! = 1, 2, or 3 (see comments above)
         depend              ! tracer dependencies (see above)

      real (kind=dbl_kind), dimension (nx_block, ny_block),   &
         intent(in) ::   &
         mflxe, mflxn   ,&! mass transport across east and north cell edges
         tarear           ! 1/tarea

      real (kind=dbl_kind), dimension (nx_block, ny_block),   &
         intent(inout) ::   &
         mm               ! mass field (mean)

      real (kind=dbl_kind), dimension (nx_block, ny_block, max_ntrace),   &
         intent(in), optional ::   &
         mtflxe, mtflxn   ! mass*tracer transport across E and N cell edges

      real (kind=dbl_kind), dimension (nx_block, ny_block, max_ntrace),   &
         intent(inout), optional ::   &
         tm               ! tracer fields

      logical (kind=log_kind), intent(inout) ::   &
         l_stop           ! if true, abort on return

      integer (kind=int_kind), intent(inout) ::   &
         istop, jstop     ! indices of grid cell where model aborts 

!
!EOP
!
      integer (kind=int_kind) ::   &
         i, j           ,&! horizontal indices
         nt, nt1, nt2     ! tracer indices

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_ntrace) ::   &
         mtold            ! old mass*tracer

      real (kind=dbl_kind) ::   &
         w1, w2           ! work variables

      integer (kind=int_kind), dimension(nx_block*ny_block) ::   &
         indxi          ,&! compressed indices in i and j directions
         indxj

      integer (kind=int_kind) ::   &
         icells         ,&! number of cells with mm > 0.
         ij               ! combined i/j horizontal index

    !-------------------------------------------------------------------
    ! Save starting values of mass*tracer
    !-------------------------------------------------------------------

      if (present(tm)) then
         do nt = 1, ntrace
            if (tracer_type(nt)==1) then ! does not depend on other tracers
               do j = jlo, jhi
               do i = ilo, ihi
                  mtold(i,j,nt) = mm(i,j) * tm(i,j,nt)
               enddo            ! i
               enddo              ! j
            elseif (tracer_type(nt)==2) then  ! depends on another tracer
               nt1 = depend(nt)
               do j = jlo, jhi
               do i = ilo, ihi
                  mtold(i,j,nt) = mm(i,j) * tm(i,j,nt1) * tm(i,j,nt)
               enddo            ! i
               enddo            ! j
            elseif (tracer_type(nt)==3) then  ! depends on two tracers
               nt1 = depend(nt)
               nt2 = depend(nt1)
               do j = jlo, jhi
               do i = ilo, ihi
                  mtold(i,j,nt) = mm(i,j)    &
                            * tm(i,j,nt2) * tm(i,j,nt1) * tm(i,j,nt)
               enddo            ! i
               enddo            ! j

 
            endif               ! depend(nt) = 0
         enddo                  ! nt
      endif                     ! present(tm)

    !-------------------------------------------------------------------
    ! Update mass field
    !-------------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi

         w1 = mflxe(i,j) - mflxe(i-1,j)   &
            + mflxn(i,j) - mflxn(i,j-1)
         mm(i,j) = mm(i,j) - w1*tarear(i,j)

         if (mm(i,j) < -puny) then    ! abort with negative value
            l_stop = .true.
            istop = i
            jstop = j
         elseif (mm(i,j) < c0) then   ! set to zero
            mm(i,j) = c0
         endif

      enddo
      enddo

      if (l_stop) then
         i = istop
         j = jstop
         w1 = mflxe(i,j) - mflxe(i-1,j)   &
            + mflxn(i,j) - mflxn(i,j-1)
         write (nu_diag,*) ' '
         write (nu_diag,*) 'New mass < 0, i, j =', i, j
         write (nu_diag,*) 'Old mass =', mm(i,j) + w1*tarear(i,j)
         write (nu_diag,*) 'New mass =', mm(i,j)
         write (nu_diag,*) 'Net transport =', -w1*tarear(i,j)
         return
      endif

    !-------------------------------------------------------------------
    ! Update tracers
    !-------------------------------------------------------------------

      if (present(tm)) then

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (mm(i,j) > c0) then ! grid cells with positive areas
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo                  ! i
         enddo                  ! j

         do nt = 1, ntrace

            do j = jlo, jhi
            do i = ilo, ihi
               tm(i,j,nt) = c0
            enddo
            enddo

            if (tracer_type(nt)==1) then ! does not depend on other tracers

               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  w1  = mtflxe(i,j,nt) - mtflxe(i-1,j,nt)   &
                      + mtflxn(i,j,nt) - mtflxn(i,j-1,nt)
                  tm(i,j,nt) = (mtold(i,j,nt) - w1*tarear(i,j))   &
                                / mm(i,j)
               enddo            ! ij


            elseif (tracer_type(nt)==2) then ! depends on another tracer
               nt1 = depend(nt)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  if (abs(tm(i,j,nt1)) > c0) then
                     w1  = mtflxe(i,j,nt) - mtflxe(i-1,j,nt)   &
                         + mtflxn(i,j,nt) - mtflxn(i,j-1,nt)
                     tm(i,j,nt) = (mtold(i,j,nt) - w1*tarear(i,j))   &
                                 / (mm(i,j) * tm(i,j,nt1))
                  endif
               enddo            ! ij

            elseif (tracer_type(nt)==3) then ! depends on two tracers
               nt1 = depend(nt)
               nt2 = depend(nt1)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  if (abs(tm(i,j,nt1)) > c0 .and.   &
                      abs(tm(i,j,nt2)) > c0) then
                     w1  = mtflxe(i,j,nt) - mtflxe(i-1,j,nt)   &
                         + mtflxn(i,j,nt) - mtflxn(i,j-1,nt)
                     tm(i,j,nt) = (mtold(i,j,nt) - w1*tarear(i,j))   &
                              / (mm(i,j) * tm(i,j,nt2) * tm(i,j,nt1))
                  endif
               enddo            ! ij
            endif               ! tracer_type
         enddo                  ! nt
      endif                     ! present(tm)

      end subroutine update_fields

!=======================================================================

      end module ice_transport_remap

!=======================================================================
