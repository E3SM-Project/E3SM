!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_remap.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!=======================================================================
!BOP
!
! !MODULE: glissade_remap - horizontal transport via incremental remapping
!
! !DESCRIPTION:
!
! Transports quantities using the second-order conservative remapping
! scheme developed by John Dukowicz and John Baumgardner (DB) and modified
! for sea ice by William Lipscomb and Elizabeth Hunke.
!
! Further modified for ice sheets by William Lipscomb.
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
! This version was created from ice_transport_remap in CICE,
!  revision 313, 6 Jan. 2011.
! The repository is here: http://oceans11.lanl.gov/svn/CICE
!
! author William H. Lipscomb, LANL
!
! !INTERFACE:
!
      module glissade_remap
!
! !USES:
!
      use glimmer_global, only: dp
      use glimmer_log
      use parallel
!
!EOP
!
      implicit none
      save
!!      private
      public :: glissade_horizontal_remap, make_remap_mask

      integer, parameter ::     &
         ngroups  = 6      ,&! number of groups of triangles that
                             ! contribute transports across each edge
         nvert = 3           ! number of vertices in a triangle

      real(dp), parameter ::   &
         puny = 1.e-11       ! small number

      !TODO - Test with bugcheck = true, but set to false for greater efficiency
      logical, parameter :: bugcheck = .true.

!TODO - Remove comments that are not relevant for CISM?
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
! The mass-like field lives in the array "mass" and the tracers fields 
! in the array "trcr". 
! In order to transport tracers correctly, the remapping routine 
! needs to know the tracers types and relationships.  This is done 
! as follows: 
! 
! Each field in the "trcr" array is assigned an index, 1:max_ntrace. 
! (Note: max_ntrace is not the same as max_ntrcr, the number of tracers 
! in the trcrn state variable array.  For remapping purposes we 
! have additional tracers hi, hs, qi and qs.) 
! For CICE with ntrcr = 1, nilyr = 4, and nslyr = 1, the 
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
! prescribed_area = T and specifying the area transported across each grid
! cell edge in the arrays edgearea_e and edgearea_n.  The departure
! regions are then tweaked, following an idea by Mats Bentsen, such
! that they have the desired area.  If prescribed_area = F, these regions
! are not tweaked, and the edgearea arrays are output variables.
!   
! Notes on the adaptation to CISM:
!
! Here we assume that all tracers are type 1, so the above arrays,
! if defined, would have the following values: 
!
!    tracer_type = (1 1 ...)
!         depend = (0 0 ...)
! has_dependents = (F F ...)
!
! But to simplify the code, these arrays have been removed throughout.
!
! Also, CISM assumes that the U grid (for velocity) is smaller than the
! T grid (for scalars).  If the T grid has dimensions (nx,ny), then the
! U grid has dimensions (nx-1,ny-1).  
!
! If nghost = 1, then there is one halo row to the south and west, but 
! there are no halo rows to the north and east.
!
! If nghost = 2, then there are two halo rows to the south and west, but 
! only one halo row to the north and east.
!
! For both the T and U grids, the local cells have dimensions (ilo:ihi,jlo:jhi),
! where ilo = 1+nghost, ihi = nx-nghost
!       jlo = 1+nghost, jhi = ny-nghost
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: glissade_horizontal_remap - incremental remapping transport scheme
!
! !INTERFACE:
!
      subroutine glissade_horizontal_remap (dt,                               &
                                            dx,                dy,            &
                                            nx_block,          ny_block,      &
                                            ntracer,           nghost,        &
                                            mmask,             icells,        &
                                            indxi,             indxj,         &
                                            uvel,              vvel,          &
                                            mass,              trcr,          &
                                            edgearea_e,        edgearea_n,    &
                                            prescribed_area_in,               &
                                            integral_order_in, dp_midpt_in)
!
! !DESCRIPTION:

! Solve the transport equations for one timestep using the incremental
! remapping scheme developed by John Dukowicz and John Baumgardner (DB)
! and modified for sea ice by William Lipscomb and Elizabeth Hunke.
!
! This scheme preserves monotonicity of mass and tracers.  That is,
! it does not produce new extrema.  It is second-order accurate in space,
! except where gradients are limited to preserve monotonicity. 
!
! This version of the remapping allows the user to specify the areal
! flux across each edge, based on an idea developed by Mats Bentsen.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!       
! !USES:
!
      use parallel
!!      use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar, horiz_bcs_stag_scalar
!
! !INPUT/OUTPUT PARAMETERS:
!
      real(dp), intent(in) ::     &
         dt           ! time step

      !TODO - Pass in dx and dy as 3D fields to allow for spatially varying
      !       cell dimensions as in POP/CICE?

      real(dp), intent(in) ::    &   
         dx, dy       ! x and y gridcell dimensions

      integer, intent(in) :: &
         nx_block   ,&! number of cells in x direction
         ny_block   ,&! number of cells in y direction
         ntracer    ,&! number of tracers to be transported
         nghost     ,&! number of ghost rows/columns 
         icells       ! number of cells with nonzero mass

      integer, intent(in), dimension(nx_block*ny_block) ::     &
         indxi, indxj     ! compressed i/j indices

      ! Note dimensions of uvel and vvel
      ! This is the CISM convention: U grid is smaller than T grid
      real(dp), intent(in), dimension(nx_block-1,ny_block-1) ::   &
         uvel       ,&! x-component of velocity (m/s)
         vvel         ! y-component of velocity (m/s)

      real(dp), intent(inout), dimension (nx_block,ny_block) ::  &
         mass       ,&! mean mass values in each grid cell
         mmask        ! = 1. if mass is present, = 0. otherwise

      real(dp), intent(inout), dimension (nx_block,ny_block,ntracer) :: &
         trcr         ! mean tracer values in each grid cell

    !-------------------------------------------------------------------
    ! If prescribed_area is true, the area of each departure region is
    !  computed in advance (e.g., by taking the divergence of the 
    !  velocity field) and passed to locate_triangles.  The departure 
    !  regions are adjusted to obtain the desired area.
    ! If false, edgearea_e and edgearea_n are computed in locate_triangles and passed out.
    !-------------------------------------------------------------------

      real(dp), dimension(nx_block,ny_block), intent(inout) ::   &
         edgearea_e     ,&! area of departure regions for east edges
         edgearea_n       ! area of departure regions for north edges

      logical, intent(in), optional ::    &
         prescribed_area_in    ! if true, edgearea_e and edgearea_n are prescribed
                               ! if false, edgearea is computed here and passed out

      integer, intent(in), optional ::    &
         integral_order_in     ! polynomial order for triangle integrals
                               ! 1 = exact for linear functions
                               ! 2 = exact for quadratic functions

      logical, intent(in), optional ::     &
         dp_midpt_in            ! if true, find departure points using
                                ! corrected midpoint velocity
!
!EOP
!
      ! local variables

      logical ::     &
         prescribed_area, dp_midpt    ! defined above

      integer :: integral_order  ! defined above

      integer ::     &
         i, j           ,&! horizontal indices
         ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      real(dp), dimension (nx_block-1,ny_block-1) ::     &
         dpx            ,&! x coordinates of departure points at cell corners
         dpy              ! y coordinates of departure points at cell corners

      real(dp), dimension(nx_block,ny_block) :: &
         mc             ,&! mass at geometric center of cell
         mx, my           ! limited derivative of mass wrt x and y

      real(dp), dimension (nx_block,ny_block,ntracer) ::     &
         tc             ,&! tracer values at geometric center of cell
         tx, ty           ! limited derivative of tracer wrt x and y

      real(dp), dimension (nx_block,ny_block) ::     &
         mflxe, mflxn     ! mass transports across E and N cell edges

      real(dp), dimension (nx_block,ny_block,ntracer) ::     &
         mtflxe, mtflxn   ! mass*tracer transports across E and N cell edges

      real(dp), dimension (nx_block,ny_block,ngroups) ::     &
         triarea          ! area of east-edge departure triangle

      real(dp), dimension (nx_block,ny_block,0:nvert,ngroups) ::  &
         xp, yp           ! x and y coordinates of special triangle points

      integer, dimension (nx_block,ny_block,ngroups) ::     &
         iflux          ,&! i index of cell contributing transport
         jflux            ! j index of cell contributing transport

      integer, dimension(ngroups) ::       &
         icellsng         ! number of cells with contribution from a given group

      integer,     &
         dimension(nx_block*ny_block,ngroups) ::     &
         indxing, indxjng ! compressed i/j indices

      logical ::     &
         l_stop           ! if true, abort the model

      character (len=5) ::   &
         edge             ! 'north' or 'east'

      real(dp), dimension(nx_block,ny_block) ::   &
          worka, workb, workc, workd

!TODO - Could save computations by passing in the following or assuming they are
!       the same for all grid cells

      real(dp), dimension (nx_block,ny_block) ::   &
         domain_mask    ,&! domain mask, = 1 wherever ice is allowed to be present
                          ! (typically = 1 everywhere)
                          ! used for gradient-limiting of mass field
         dxt            ,&! T-cell width (m)
         dyt            ,&! T-cell height (m)
         dxu            ,&! U-cell width (m)
         dyu            ,&! U-cell height (m)
         htn            ,&! length of north cell edge (m)
         hte              ! length of east cell edge (m)

      real(dp) ::   &
         tarear           ! reciprocal grid cell area

      real(dp), dimension (nx_block,ny_block) ::   &
         xav, yav          ,&! gridcell avg values of x, y
         xxav, xyav, yyav    ! gridcell avg values of xx, xy, yy

      character(len=100) :: message

    !------------------------------------------------------------------- 
    ! Initialize various grid quantities and code options
    !------------------------------------------------------------------- 

      ! Assume that ice can exist everywhere on the domain
      ! May need to pass this in as an argument if parts of the domain are masked out,
      !  e.g. Ellesmere Island for a Greenland ice sheet simulation.

      domain_mask(:,:) = 1.d0

      ! Assume gridcells are rectangular, in which case dxt = dxu = htn
      !  and dyt = dyu = hte.
      !TODO - Pass in dx and dy as 3D fields to allow for spatially varying
      !       cell dimensions as in POP/CICE?

      dxt(:,:) = dx
      dxu(:,:) = dx
      htn(:,:) = dx

      dyt(:,:) = dy
      dyu(:,:) = dy
      hte(:,:) = dy

      if (dx*dy > 0.d0) then
         tarear = 1.d0 / (dx*dy)
      else
         tarear = 0.d0
      endif

      xav(:,:) = 0.d0
      yav(:,:) = 0.d0
      xxav(:,:) = 1.d0 / 12.d0  ! These are the scaled values, valid if dxt = dyt = 1
      yyav(:,:) = 1.d0 / 12.d0
!!      xxav(:,:) = dxt(:,:)**2 / 12.d0  ! These would be used if dimensional values 
!!      yyav(:,:) = dyt(:,:)**2 / 12.d0  ! of dxt and dyt were passed to construct_fields
      xyav(:,:) = 0.d0

      l_stop = .false.

      if (present(integral_order_in)) then
         integral_order = integral_order_in
      else
         integral_order = 2  ! exact for integrating quadratic functions
      endif

      if (present(dp_midpt_in)) then
         dp_midpt = dp_midpt_in
      else
!WHL - Set to true for increased accuracy
!    - Set to false for closer agreement with the old remapping code
!!!         dp_midpt = .true.
!pw++
         dp_midpt = .false.
!pw--
      endif

      if (present(prescribed_area_in)) then
         prescribed_area = prescribed_area_in
      else
         prescribed_area = .false.
      endif

      ! These arrays are passed to construct_fields in lieu of the dimensional
      ! values of dxt, dyt, htn, and hte.
      ! 
      worka(:,:) = 1.d0
      workb(:,:) = 1.d0
      workc(:,:) = 1.d0
      workd(:,:) = 1.d0

      ! Compute lower and upper indices for locally owned cells
      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

!      print*, 'ilo, ihi =', ilo, ihi
!      print*, 'jlo, jhi =', jlo, jhi

    !-------------------------------------------------------------------
    ! Construct linear fields, limiting gradients to preserve monotonicity.
    ! Note: Pass in unit arrays instead of true distances hte, htn, etc.
    !       The resulting gradients are in scaled coordinates.
    !-------------------------------------------------------------------

      call construct_fields(nx_block,            ny_block,           &
                            ilo, ihi,            jlo, jhi,           &
                            nghost,              ntracer,            &
                            icells,                                  &
                            indxi    (:),        indxj(:),           &
!                            htn      (:,:),      hte    (:,:),      &
                            worka    (:,:),      workb  (:,:),       &
                            domain_mask(:,:),    xav    (:,:),       &
                            yav      (:,:),      xxav   (:,:),       &
                            xyav     (:,:),      yyav   (:,:),       &
!                            dxt      (:,:),      dyt  (:,:),        &
                            workc    (:,:),      workd  (:,:),       &
                            mass  (:,:),         mc  (:,:),          &
                            mx    (:,:),         my  (:,:),          &
                            mmask (:,:),                             &
                            trcr  (:,:,:),       tc  (:,:,:),        &
                            tx    (:,:,:),       ty  (:,:,:))

    !-------------------------------------------------------------------
    ! Given velocity field at cell corners, compute departure points
    ! of trajectories.
    !-------------------------------------------------------------------

      call departure_points(nx_block,         ny_block,          &
                            ilo, ihi,         jlo, jhi,          &
                            nghost,           dt,                &
                            uvel  (:,:),      vvel(:,:),         &
                            dxu   (:,:),      dyu (:,:),         &
                            htn   (:,:),      hte (:,:),         &
                            dpx   (:,:),      dpy (:,:),         &
                            dp_midpt,         l_stop)

      if (l_stop) then
         write(message,*) 'Aborting (task = ',this_rank,')'
         call write_log(message,GM_FATAL)
      endif

    !-------------------------------------------------------------------
    ! Ghost cell updates
    ! If nghost >= 2, these calls are not needed
    !-------------------------------------------------------------------

      if (nghost==1) then

         ! mass field
         call parallel_halo(mc)
!!         call horiz_bcs_unstag_scalar(mc)
         call parallel_halo(mx)
!!         call horiz_bcs_unstag_scalar(mx)
         call parallel_halo(my)
!!         call horiz_bcs_unstag_scalar(my)

         ! tracer fields
         call parallel_halo(tc)
!!         call horiz_bcs_unstag_scalar(tc)
         call parallel_halo(tx)
!!         call horiz_bcs_unstag_scalar(tx)
         call parallel_halo(ty)
!!         call horiz_bcs_unstag_scalar(ty)

         ! departure points
         call parallel_halo(dpx)
!!         call horiz_bcs_stag_scalar(dpx)
         call parallel_halo(dpy)
!!         call horiz_bcs_stag_scalar(dpy)

      endif  ! nghost

    !-------------------------------------------------------------------
    ! Transports for east cell edges.
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Compute areas and vertices of departure triangles.
    !-------------------------------------------------------------------

      edge = 'east'
      call locate_triangles(nx_block,          ny_block,           &
                            ilo, ihi,          jlo, jhi,           &
                            edge,              icellsng (:),       &
                            indxing(:,:),      indxjng(:,:),       &
                            dpx  (:,:),        dpy (:,:),          &
                            dxu  (:,:),        dyu (:,:),          &
                            xp(:,:,:,:),       yp(:,:,:,:),        &
                            iflux(:,:,:),      jflux(:,:,:),       &
                            triarea(:,:,:),                        &
                            prescribed_area,   edgearea_e(:,:))     

    !-------------------------------------------------------------------
    ! Given triangle vertices, compute coordinates of triangle points
    !  needed for transport integrals.
    !-------------------------------------------------------------------

      call triangle_coordinates (nx_block,          ny_block,      &
                                 icellsng (:),                     &
                                 indxing(:,:),      indxjng(:,:),  &
                                 xp,                yp,            &
                                 integral_order)

    !-------------------------------------------------------------------
    ! Compute the transport across east cell edges by summing contributions
    ! from each triangle.
    !-------------------------------------------------------------------

      call transport_integrals(nx_block,          ny_block,           &
                               ntracer,           icellsng (:),       &
                               indxing(:,:),      indxjng(:,:),       &
                               triarea(:,:,:),    integral_order,     &
                               iflux(:,:,:),      jflux(:,:,:),       &
                               xp(:,:,:,:),       yp(:,:,:,:),       &
                               mc(:,:),           mx   (:,:),         &
                               my(:,:),           mflxe(:,:),         &
                               tc(:,:,:),         tx   (:,:,:),       &
                               ty(:,:,:),         mtflxe(:,:,:))

    !-------------------------------------------------------------------
    ! Repeat for north edges
    !-------------------------------------------------------------------

      edge = 'north'
      call locate_triangles(nx_block,          ny_block,           &
                            ilo, ihi,          jlo, jhi,           &
                            edge,              icellsng (:),       &
                            indxing(:,:),      indxjng(:,:),       &
                            dpx  (:,:),        dpy (:,:),          &
                            dxu  (:,:),        dyu (:,:),          &
                            xp(:,:,:,:),       yp(:,:,:,:),        &
                            iflux(:,:,:),      jflux(:,:,:),       &
                            triarea(:,:,:),                        &
                            prescribed_area,   edgearea_n(:,:))     

      call triangle_coordinates (nx_block,          ny_block,      &
                                 icellsng (:),                     &
                                 indxing(:,:),      indxjng(:,:),  &
                                 xp,                yp,            &
                                 integral_order)

      call transport_integrals(nx_block,           ny_block,          &
                               ntracer,            icellsng (:),      &
                               indxing(:,:),       indxjng(:,:),      &
                               triarea(:,:,:),     integral_order,    &
                               iflux(:,:,:),       jflux(:,:,:),      &
                               xp(:,:,:,:),        yp(:,:,:,:),       &
                               mc(:,:),            mx    (:,:),       &
                               my(:,:),            mflxn (:,:),       &
                               tc(:,:,:),          tx    (:,:,:),     &
                               ty(:,:,:),          mtflxn(:,:,:))

    !-------------------------------------------------------------------
    ! Update the ice area and tracers.
    !-------------------------------------------------------------------

      call update_fields (nx_block,           ny_block,          &
                          ilo, ihi,           jlo, jhi,          &
                          ntracer,                               &
                          tarear,             l_stop,            &
                          mflxe (:,:),        mflxn (:,:),       &
                          mass  (:,:),                           &     
                          mtflxe(:,:,:),      mtflxn(:,:,:),     &
                          trcr  (:,:,:) )

      if (l_stop) then
         write(message,*) 'Aborting (task = ',this_rank,')'
         call write_log(message,GM_FATAL)
      endif

      end subroutine glissade_horizontal_remap

!=======================================================================
!
!BOP
!
! !IROUTINE: make_remap_mask - make ice mask
!
! !INTERFACE:
!
      subroutine make_remap_mask (nx_block, ny_block,           &
                                  ilo, ihi, jlo, jhi,           &
                                  nghost,   icells,             &
                                  indxi,    indxj,              &
                                  mass,     mmask)
!
! !DESCRIPTION:
!
! Make ice mask; identify cells where ice is present.
!
! If a gridcell is massless (mass < puny), then the values of tracers
!  in that grid cell are assumed to have no physical meaning.
!WHL - Changed this condition from 'mass < puny' to 'mass < 0.d0'
!      to preserve monotonicity in grid cells with very small thickness
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer, intent(in) ::     &
           nx_block, ny_block  ,&! block dimensions
           ilo,ihi,jlo,jhi     ,&! beginning and end of physical domain
           nghost                ! number of ghost cells

      integer, intent(out) ::     &
           icells         ! number of cells with ice

      integer, dimension(nx_block*ny_block), intent(out) ::     &
           indxi        ,&! compressed i/j indices
           indxj

      real(dp), dimension (nx_block,ny_block),            &
           intent(in) ::     &
           mass          ! mean ice thickness in each grid cell

      real(dp), dimension (nx_block,ny_block),            &
           intent(out) ::     &
           mmask         ! = 1. if ice is present, else = 0.
!
!EOP
!
      integer ::     &
           i, j, ij      ! indices

    !-------------------------------------------------------------------
    ! ice mask
    !-------------------------------------------------------------------

!WHL - Changed this condition from 'mass(i,j) < puny' to 'mass(i,j) < 0.d0'
!      to preserve monotonicity in grid cells with very small thickness
      do j = 1, ny_block
      do i = 1, nx_block
!!         if (mass(i,j) > puny) then
         if (mass(i,j) > 0.d0) then
            mmask(i,j) = 1.d0
         else
            mmask(i,j) = 0.d0
         endif
      enddo
      enddo

    !-------------------------------------------------------------------
    ! Tag grid cells where ice is present
    ! For nghost = 1, exclude ghost cells
    ! For nghost = 2, include one layer of ghost cells
    !-------------------------------------------------------------------

      icells = 0
      do ij = 1, nx_block*ny_block
        indxi(ij) = 0
        indxj(ij) = 0
      enddo

      do j = jlo-nghost+1, jhi+nghost-1
      do i = ilo-nghost+1, ihi+nghost-1
!WHL - Changed this condition from 'mass(i,j) > puny' to 'mass(i,j) > 0.d0'
!      to preserve monotonicity in grid cells with very small thickness
!!         if (mass(i,j) > puny) then
         if (mass(i,j) > 0.d0) then
            icells = icells + 1
            ij = icells
            indxi(ij) = i
            indxj(ij) = j
         endif
      enddo
      enddo

      end subroutine make_remap_mask

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
                                   nghost,         ntracer,    &
                                   icells,                     &
                                   indxi,          indxj,      &
                                   htn,            hte,        &
                                   hm,             xav,        &
                                   yav,            xxav,       &
                                   xyav,           yyav,       &
                                   dxt,            dyt,        &
                                   mass,           mc,         &
                                   mx,             my,         &
                                   mmask,                      &
                                   trcr,           tc,         &
                                   tx,             ty)
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
      integer, intent(in) ::   &
         nx_block, ny_block  ,&! block dimensions
         ilo,ihi,jlo,jhi     ,&! beginning and end of physical domain
         nghost              ,&! number of ghost cell layers
         ntracer             ,&! number of tracers in use
         icells                ! number of cells with mass

      integer, dimension(nx_block*ny_block), intent(in) :: &
         indxi          ,&! compressed i/j indices
         indxj

      real(dp), dimension (nx_block,ny_block), intent(in) ::   &
         hm             ,&! domain mask
         htn            ,&! length of northern edge of T-cell (m)
         hte            ,&! length of eastern edge of T-cell (m)
         xav,  yav      ,&! mean T-cell values of x, y
         xxav, xyav, yyav ,&! mean T-cell values of xx, xy, yy
         dxt            ,&! grid cell width (m)
         dyt              ! grid cell height (m)

      real(dp), dimension (nx_block,ny_block), intent(in) ::   &
         mass          ,&! mean value of mass field
         mmask           ! = 1. if ice is present, = 0. otherwise

      real(dp), dimension (nx_block,ny_block), intent(out) ::   &
         mc             ,&! mass value at geometric center of cell
         mx, my           ! limited derivative of mass wrt x and y

      real(dp), dimension (nx_block,ny_block,ntracer),  &
         intent(in), optional ::   &
         trcr             ! mean tracer

      real(dp), dimension (nx_block,ny_block,ntracer),  &
         intent(out), optional ::   &
         tc             ,&! tracer at geometric center of cell
         tx, ty           ! limited derivative of tracer wrt x and y
!
!EOP
!
      integer ::   &
         i, j           ,&! horizontal indices
         nt, nt1        ,&! tracer indices
         ij               ! combined i/j horizontal index

      real(dp), dimension (nx_block,ny_block) ::    &
         mxav           ,&! x coordinate of center of mass
         myav             ! y coordinate of center of mass

      real(dp), dimension (nx_block,ny_block,ntracer) ::  &
         mtxav          ,&! x coordinate of center of mass*tracer
         mtyav            ! y coordinate of center of mass*tracer

      real(dp) ::   &
         w1, w2, w3, w4, w5, w6, w7   ! work variables

    !-------------------------------------------------------------------
    ! Compute field values at the geometric center of each grid cell,
    ! and compute limited gradients in the x and y directions.
    !
    ! For second order accuracy, each state variable is approximated as
    ! a field varying linearly over x and y within each cell.  For each
    ! category, the integrated value of m(x,y) over the cell must
    ! equal mass(i,j,n)*tarea(i,j), where tarea(i,j) is the cell area.
    ! Similarly, the integrated value of m(x,y)*t(x,y) must equal
    ! the total mass*tracer, mass(i,j,n)*trcr(i,j,n)*tarea(i,j).
    !
    ! These integral conditions are satisfied for linear fields if we
    ! stipulate the following:
    ! (1) The mean mass is equal to the mass at the cell centroid.
    ! (2) The mean value trcr1 of type 1 tracers is equal to the value
    !     at the center of mass.
    ! (3) The mean value trcr2 of type 2 tracers is equal to the value
    !     at the center of mass*trcr1, where trcr2 depends on trcr1.
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
         mc(i,j)  = 0.d0
         mx(i,j)  = 0.d0
         my(i,j)  = 0.d0
         mxav(i,j) = 0.d0
         myav(i,j) = 0.d0
      enddo
      enddo

      if (present(trcr)) then
         do nt = 1, ntracer
           do j = 1, ny_block
            do i = 1, nx_block
               tc(i,j,nt) = 0.d0
               tx(i,j,nt) = 0.d0
               ty(i,j,nt) = 0.d0
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
                             mass,     hm,         &
                             xav,      yav,        &
                             htn,      hte,        &
                             dxt,      dyt,        &
                             mx,       my)

      do ij = 1,icells   ! ice is present
         i = indxi(ij)
         j = indxj(ij)

         ! mass field at geometric center
         mc(i,j) = mass(i,j) - xav(i,j)*mx(i,j)   &
                             - yav(i,j)*my(i,j)

      enddo                     ! ij

      ! tracers

      if (present(trcr)) then

       do ij = 1,icells       ! cells with mass
          i = indxi(ij)
          j = indxj(ij)

         ! center of mass (mxav,myav) for each cell
          mxav(i,j) = (mx(i,j)*xxav(i,j)    &
                     + my(i,j)*xyav(i,j)    &
                     + mc(i,j)*xav (i,j)) / mass(i,j)
          myav(i,j) = (mx(i,j)*xyav(i,j)    &
                     + my(i,j)*yyav(i,j)    &
                     + mc(i,j)*yav (i,j)) / mass(i,j)
       enddo

       do nt = 1, ntracer

         call limited_gradient(nx_block,     ny_block,  &
                                ilo, ihi,     jlo, jhi,  &
                                nghost,                  &
                                trcr(:,:,nt), mmask,     &
                                mxav,         myav,      &
                                htn,          hte,       &
                                dxt,          dyt,       &
                                tx(:,:,nt),   ty(:,:,nt)) 

          do ij = 1, icells      ! mass is present
             i = indxi(ij)
             j = indxj(ij)

             ! tracer value at geometric center
             tc(i,j,nt) = trcr(i,j,nt) - tx(i,j,nt)*mxav(i,j)   &
                                       - ty(i,j,nt)*myav(i,j)
          enddo            ! ij

       enddo                    ! ntracer

     endif                     ! present (trcr)

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
                                   htn,      hte,        &
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
      integer, intent(in) ::   &
          nx_block, ny_block,&! block dimensions
          ilo,ihi,jlo,jhi ,&! beginning and end of physical domain
          nghost              ! number of ghost cell layers

      real(dp), dimension (nx_block,ny_block),   &
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
          htn    ,&! length of northern edge of T-cell (m)
          hte      ! length of eastern edge of T-cell (m)

      real(dp), dimension (nx_block,ny_block),   &
          intent(out) ::   &
          gx     ,&! limited x-direction gradient
          gy       ! limited y-direction gradient
!
!EOP
!
      integer ::   &
          i, j, ij        ,&! standard indices
          icells            ! number of cells to limit

      integer, dimension(nx_block*ny_block) ::   &
          indxi, indxj   ! combined i/j horizontal indices

      real(dp) ::   &
          phi_nw, phi_n, phi_ne ,&! values of phi in 8 neighbor cells
          phi_w,         phi_e  ,&
          phi_sw, phi_s, phi_se ,&
          qmn, qmx     ,&! min and max value of phi within grid cell
          pmn, pmx     ,&! min and max value of phi among neighbor cells
          w1, w2, w3, w4 ! work variables

      real(dp) ::   &
          gxtmp, gytmp   ! temporary term for x- and y- limited gradient

      gx(:,:) = 0.d0
      gy(:,:) = 0.d0

      ! For nghost = 1, loop over physical cells and update ghost cells later
      ! For nghost = 2, loop over a layer of ghost cells and skip the update

      icells = 0
      do j = jlo-nghost+1, jhi+nghost-1
      do i = ilo-nghost+1, ihi+nghost-1
         if (phimask(i,j) > puny) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif                  ! phimask > puny
      enddo
      enddo

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! Store values of phi in the 8 neighbor cells.
         ! Note: phimask = 1. or 0.  If phimask = 1., use the true value;
         !  if phimask = 0., use the home cell value so that non-physical
         !  values of phi do not contribute to the gradient.
         phi_nw = phimask(i-1,j+1) * phi(i-1,j+1)   &
            + (1.d0-phimask(i-1,j+1))* phi(i,j)
         phi_n  = phimask(i,j+1)   * phi(i,j+1)   &
            + (1.d0-phimask(i,j+1))  * phi(i,j)
         phi_ne = phimask(i+1,j+1) * phi(i+1,j+1)   &
            + (1.d0-phimask(i+1,j+1))* phi(i,j)
         phi_w  = phimask(i-1,j)   * phi(i-1,j)   &
            + (1.d0-phimask(i-1,j))  * phi(i,j)
         phi_e  = phimask(i+1,j)   * phi(i+1,j)   &
            + (1.d0-phimask(i+1,j))  * phi(i,j)
         phi_sw = phimask(i-1,j-1) * phi(i-1,j-1)   &
            + (1.d0-phimask(i-1,j-1))* phi(i,j)
         phi_s  = phimask(i,j-1)   * phi(i,j-1)   &
            + (1.d0-phimask(i,j-1))  * phi(i,j)
         phi_se = phimask(i+1,j-1) * phi(i+1,j-1)   &
            + (1.d0-phimask(i+1,j-1))* phi(i,j)

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

         w1  =  (0.5d0*htn(i,j)   - cnx(i,j)) * gxtmp   &
              + (0.5d0*hte(i,j)   - cny(i,j)) * gytmp
         w2  =  (0.5d0*htn(i,j-1) - cnx(i,j)) * gxtmp   &
              - (0.5d0*hte(i,j)   + cny(i,j)) * gytmp
         w3  = -(0.5d0*htn(i,j-1) + cnx(i,j)) * gxtmp   &
              - (0.5d0*hte(i-1,j) + cny(i,j)) * gytmp
         w4  =  (0.5d0*hte(i-1,j) - cny(i,j)) * gytmp   &
              - (0.5d0*htn(i,j)   + cnx(i,j)) * gxtmp

         qmn = min (w1, w2, w3, w4)
         qmx = max (w1, w2, w3, w4)

         ! the limiting coefficient
         if (abs(qmn) > 0.d0) then ! 'abs(qmn) > puny' not sufficient
            w1 = max(0.d0, pmn/qmn)
         else
            w1 = 1.d0
         endif

         if (abs(qmx) > 0.d0) then
            w2 = max(0.d0, pmx/qmx)
         else
            w2 = 1.d0
         endif

         w1 = min(1.d0, w1, w2)

         ! Limit the gradient components
         gx(i,j) = w1 * gxtmp
         gy(i,j) = w1 * gytmp

      enddo                     ! ij

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
                                   htn,        hte,     &
                                   dpx,        dpy,     &
                                   dp_midpt,   l_stop)
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
      integer, intent(in) ::   &
         nx_block, ny_block,&! block dimensions
         ilo,ihi,jlo,jhi,   &! beginning and end of physical domain
         nghost              ! number of ghost cell layers

      real(dp), intent(in) ::   &
         dt               ! time step (s)

      real(dp), dimension (nx_block-1,ny_block-1), intent(in) ::   &
         uvel           ,&! x-component of velocity (m/s)
         vvel             ! y-component of velocity (m/s)

      real(dp), dimension (nx_block-1,ny_block-1), intent(out) ::   &
         dpx            ,&! coordinates of departure points (m)
         dpy              ! coordinates of departure points (m)

      real(dp), dimension (nx_block,ny_block), intent(in) ::   &
         dxu            ,&! E-W dimensions of U-cell (m)
         dyu            ,&! N-S dimensions of U-cell (m)
         htn            ,&! length of north face of T-cell (m) 
         hte              ! length of east face of T-cell (m) 

      logical, intent(in) ::   &
         dp_midpt            ! if true, find departure points using
                             ! corrected midpoint velocity

      logical, intent(inout) ::   &
         l_stop       ! if true, abort on return
!
!EOP
!
      integer ::   &
         i, j, i2, j2     ! horizontal indices

      real(dp) ::                  &
         mpx,  mpy      ,&! coordinates of midpoint of back trajectory,
                          ! relative to cell corner
         mpxt, mpyt     ,&! midpoint coordinates relative to cell center
         ump,  vmp        ! corrected velocity at midpoint

      integer ::   &
         istop, jstop     ! indices of grid cell where model aborts 

      character(len=100) :: message

    !-------------------------------------------------------------------
    ! Estimate departure points.
    ! This estimate is 1st-order accurate in time; improve accuracy by
    !  using midpoint approximation (to add later).
    ! For nghost = 1, loop over physical cells and update ghost cells later.
    ! For nghost = 2, loop over a layer of ghost cells and skip update.
    !-------------------------------------------------------------------

      dpx(:,:) = 0.d0
      dpy(:,:) = 0.d0

      ! Note: If nghost = 1, then this loop will include all vertices of all locally owned cells,
      !         including halo values along the west and south edges of the domain.
      !       If nghost = 2, then this loop includes an additional layer of cells around the domain,
      !         as needed if we are using the midpoint correction method.

      do j = jlo-nghost, jhi+nghost-1
      do i = ilo-nghost, ihi+nghost-1

         dpx(i,j) = -dt*uvel(i,j)
         dpy(i,j) = -dt*vvel(i,j)

         ! Check for values out of bounds (more than one grid cell away)
         if (dpx(i,j) < -htn(i,j) .or. dpx(i,j) > htn(i+1,j) .or.   &
             dpy(i,j) < -hte(i,j) .or. dpy(i,j) > hte(i,j+1)) then

!WHL - debug
!             print*, ' '
!             print*, 'dt =', dt
!             print*, 'i, j =', i, j
!             print*, 'dpx, dpy =', dpx(i,j), dpy(i,j)
!             print*, 'hte, htn =', hte(i,j), htn(i,j)
!             print*, 'bad departure points'

            l_stop = .true.
            istop = i
            jstop = j
         endif

      enddo
      enddo

!TODO - Write error message cleanly to the log file.
!       I think this will require broadcasting istop and jstop to main_rank.
!       For now, just print an error message locally.

      if (l_stop) then
         i = istop
         j = jstop
!         write (message,*) 'Process:',this_rank
!         call write_log(message)
!         write (message,*) 'Remap, departure points out of bounds:, i, j =', i, j
!         call write_log(message)
!         write (message,*) 'dpx, dpy =', dpx(i,j), dpy(i,j)
!         call write_log(message)
!         write (message,*) 'uvel, vvel =', uvel(i,j), vvel(i,j)
!         call write_log(message)
!         write (message,*) 'htn(i,j), htn(i+1,j) =', htn(i,j), htn(i+1,j)
!         call write_log(message)
!         write (message,*) 'hte(i,j), hte(i,j+1) =', hte(i,j), hte(i,j+1)
!         call write_log(message)
         write (6,*) 'Process:', this_rank
         write (6,*) 'Remap, departure points out of bounds:, i, j =', i, j
         write (6,*) 'dpx, dpy =', dpx(i,j), dpy(i,j)
         write (6,*) 'uvel, vvel =', uvel(i,j), vvel(i,j)
         write (6,*) 'htn(i,j), htn(i+1,j) =', htn(i,j), htn(i+1,j)
         write (6,*) 'hte(i,j), hte(i,j+1) =', hte(i,j), hte(i,j+1)
         return
      endif

!Note: Need nghost >= 2 to do this correction, which requires velocities
!      for vertices with indices (ilo-2) and (jlo-2).
 
      if (dp_midpt .and. nghost>= 2) then   ! find dep pts using corrected midpt velocity 

       do j = jlo-1, jhi
       do i = ilo-1, ihi

         if (uvel(i,j)/=0.d0 .or. vvel(i,j)/=0.d0) then
 
    !-------------------------------------------------------------------
    ! Scale departure points to coordinate system in which grid cells
    ! have sides of unit length.
    !-------------------------------------------------------------------

            dpx(i,j) = dpx(i,j) / dxu(i,j)
            dpy(i,j) = dpy(i,j) / dyu(i,j)

    !-------------------------------------------------------------------
    ! Estimate midpoint of backward trajectory relative to corner (i,j).
    !-------------------------------------------------------------------

            mpx = 0.5d0 * dpx(i,j)
            mpy = 0.5d0 * dpy(i,j)
 
    !-------------------------------------------------------------------
    ! Determine the indices (i2,j2) of the cell where the trajectory lies.
    ! Compute the coordinates of the midpoint of the backward trajectory
    !  relative to the cell center in a stretched coordinate system
    !  with vertices at (1/2, 1/2), (1/2, -1/2), etc.
    !-------------------------------------------------------------------

            if (mpx >= 0.d0 .and. mpy >= 0.d0) then    ! cell (i+1,j+1)
               i2 = i+1
               j2 = j+1
               mpxt = mpx - 0.5d0
               mpyt = mpy - 0.5d0
            elseif (mpx < 0.d0 .and. mpy < 0.d0) then  ! cell (i,j)
               i2 = i
               j2 = j
               mpxt = mpx + 0.5d0
               mpyt = mpy + 0.5d0
            elseif (mpx >= 0.d0 .and. mpy < 0.d0) then ! cell (i+1,j)
               i2 = i+1
               j2 = j
               mpxt = mpx - 0.5d0
               mpyt = mpy + 0.5d0
            elseif (mpx < 0.d0 .and. mpy >= 0.d0) then ! cell (i,j+1)
               i2 = i
               j2 = j+1
               mpxt = mpx + 0.5d0
               mpyt = mpy - 0.5d0
            endif
            
    !-------------------------------------------------------------------
    ! Using a bilinear approximation, estimate the velocity at the
    ! trajectory midpoint in the (i2,j2) reference frame.
    !-------------------------------------------------------------------
 
            ump = uvel(i2-1,j2-1)*(mpxt-0.5d0)*(mpyt-0.5d0)     &
                - uvel(i2,  j2-1)*(mpxt+0.5d0)*(mpyt-0.5d0)     &
                + uvel(i2,  j2  )*(mpxt+0.5d0)*(mpyt+0.5d0)     &  
                - uvel(i2-1,j2  )*(mpxt-0.5d0)*(mpyt+0.5d0)
 
            vmp = vvel(i2-1,j2-1)*(mpxt-0.5d0)*(mpyt-0.5d0)     &
                - vvel(i2,  j2-1)*(mpxt+0.5d0)*(mpyt-0.5d0)     &
                + vvel(i2,  j2  )*(mpxt+0.5d0)*(mpyt+0.5d0)     &
                - vvel(i2-1,j2  )*(mpxt-0.5d0)*(mpyt+0.5d0)
 
    !-------------------------------------------------------------------
    ! Use the midpoint velocity to estimate the coordinates of the
    !  departure point relative to corner (i,j).
    !-------------------------------------------------------------------
 
            dpx(i,j) = -dt * ump
            dpy(i,j) = -dt * vmp
 
         endif               ! nonzero velocity

       enddo                 ! i
       enddo                 ! j
 
      endif                  ! dp_midpt

      end subroutine departure_points

!=======================================================================
!
!BOP
!
! !IROUTINE: locate_triangles - triangle info for cell edges
!
! !INTERFACE:
!
      subroutine locate_triangles (nx_block,        ny_block,   &
                                   ilo, ihi,        jlo, jhi,   &
                                   edge,            icells,     &
                                   indxi,           indxj,      &
                                   dpx,             dpy,        &
                                   dxu,             dyu,        &
                                   xp,              yp,         &
                                   iflux,           jflux,      &
                                   triarea,                  &
                                   prescribed_area, edgearea)
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
      integer, intent(in) ::   &
         nx_block, ny_block,&! block dimensions
         ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      character (len=5), intent(in) ::   &
         edge             ! 'north' or 'east'

      real(dp), dimension(nx_block-1,ny_block-1), intent(in) ::  &
         dpx            ,&! x coordinates of departure points at cell corners
         dpy              ! y coordinates of departure points at cell corners

      real(dp), dimension(nx_block,ny_block), intent(in) ::  &
         dxu            ,&! E-W dimension of U-cell (m)
         dyu              ! N-S dimension of U-cell (m)

      real(dp), dimension (nx_block,ny_block,0:nvert,ngroups),   &
         intent(out) ::   &
         xp, yp           ! coordinates of triangle vertices

      real(dp), dimension (nx_block,ny_block,ngroups),   &
           intent(out) ::   &
         triarea          ! area of departure triangle

      integer, dimension (nx_block,ny_block,ngroups),    &
         intent(out) ::   &
         iflux          ,&! i index of cell contributing transport
         jflux            ! j index of cell contributing transport

      integer, dimension (ngroups), intent(out) ::   &
         icells           ! number of cells where triarea > puny

      integer, dimension (nx_block*ny_block,ngroups), &
         intent(out) ::                                               &
         indxi          ,&! compressed index in i-direction
         indxj            ! compressed index in j-direction

      logical, intent(in) ::   &
         prescribed_area  ! if true, the area of each departure region is
                          !  passed in as edgearea
                          ! if false, edgearea if determined internally
                          !  and is passed out
                          
      real(dp), dimension(nx_block,ny_block), intent(inout) ::   &
         edgearea         ! area of departure region for each edge
                          ! edgearea > 0 for eastward/northward flow
!
!EOP
!
      integer ::   &
         i, j, ij, ic   ,&! horizontal indices
         ib, ie, jb, je ,&! limits for loops over edges
         ng, nv         ,&! triangle indices
         ishift, jshift   ! differences between neighbor cells

      integer ::   &
         icellsd          ! number of cells where departure area > 0.

      integer, dimension (nx_block*ny_block) ::  &
         indxid         ,&! compressed index in i-direction
         indxjd           ! compressed index in j-direction

      real(dp), dimension(nx_block,ny_block) ::   &
         dx, dy         ,&! scaled departure points
         areafac_c      ,&! area scale factor at center of edge
         areafac_l      ,&! area scale factor at left corner
         areafac_r        ! area scale factor at right corner

      real(dp) ::   &
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
                          ! shifted if prescribed_area = T
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

      real(dp), dimension (nx_block,ny_block,ngroups) ::   &
         areafact         ! = 1 for positive flux, -1 for negative

      real(dp), dimension(nx_block,ny_block) ::   &
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
    ! (4) BC1a, BC1b, TC1a, TC1b
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
    ! and the transport is across the edge between cells TC and BC.
    !
    ! Departure points are scaled to a local coordinate system
    !  whose origin is at the midpoint of the edge.
    ! In this coordinate system, the lefthand corner CL = (-0.5,0)
    !  and the righthand corner CR = (0.5, 0).
    !-------------------------------------------------------------------
  
    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      dx(:,:) = 0.d0
      dy(:,:) = 0.d0
      areafac_c(:,:) = 0.d0
      areafac_l(:,:) = 0.d0
      areafac_r(:,:) = 0.d0
      do ng = 1, ngroups
         do j = 1, ny_block
         do i = 1, nx_block
            triarea (i,j,ng) = 0.d0
            areafact(i,j,ng) = 0.d0
            iflux   (i,j,ng) = i
            jflux   (i,j,ng) = j
         enddo
         enddo
         do nv = 0, nvert
            do j = 1, ny_block
            do i = 1, nx_block
               xp(i,j,nv,ng) = 0.d0
               yp(i,j,nv,ng) = 0.d0
            enddo
            enddo
         enddo
      enddo

      if (trim(edge) == 'north') then

         ! loop size
         ! Note: The loop is over all north edges that border one or more locally owned grid
         !  cells. This includes the north edges of cells with index (jlo-1), which are the south
         !  edges of cells with index (jlo).

         ib = ilo
         ie = ihi 
         jb = jlo - 1            ! lowest j index is a ghost cell
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
            areafac_c(i,j) = 0.5d0*(areafac_l(i,j) + areafac_r(i,j))
         enddo
         enddo

      else                      ! east edge

         ! loop size
         ! Note: The loop is over all east edges that border one or more locally owned grid
         !  cells. This includes the east edges of cells with index (ilo-1), which are the west
         !  edges of cells with index (ilo).

         ib = ilo - 1            ! lowest i index is a ghost cell
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
            areafac_c(i,j) = 0.5d0 * (areafac_l(i,j) + areafac_r(i,j))
         enddo
         enddo

      endif

    !-------------------------------------------------------------------
    ! Compute mask for edges with nonzero departure areas
    !-------------------------------------------------------------------

      if (prescribed_area) then
         icellsd = 0
         do j = jb, je
         do i = ib, ie
            if (edgearea(i,j) /= 0.d0) then
               icellsd = icellsd + 1
               indxid(icellsd) = i
               indxjd(icellsd) = j
            endif
         enddo
         enddo
      else
         icellsd = 0
         if (trim(edge) == 'north') then
            do j = jb, je   ! jb = jlo - 1
            do i = ib, ie   ! ib = ilo
               if (dpx(i-1,j)/=0.d0 .or. dpy(i-1,j)/=0.d0   &
                                  .or.                  &
                     dpx(i,j)/=0.d0 .or.   dpy(i,j)/=0.d0) then
                  icellsd = icellsd + 1
                  indxid(icellsd) = i
                  indxjd(icellsd) = j
               endif
            enddo
            enddo
         else       ! east edge
            do j = jb, je   ! jb = jlo
            do i = ib, ie   ! ib = ilo - 1
               if (dpx(i,j-1)/=0.d0 .or. dpy(i,j-1)/=0.d0   &
                                  .or.                  &
                     dpx(i,j)/=0.d0 .or.   dpy(i,j)/=0.d0) then
                  icellsd = icellsd + 1
                  indxid(icellsd) = i
                  indxjd(icellsd) = j
               endif
            enddo
            enddo
         endif       ! edge = north/east
      endif          ! prescribed_area

    !-------------------------------------------------------------------
    ! Scale the departure points.
    ! Note: This loop must include all vertices of all edges for which
    !       fluxes are computed.
    !-------------------------------------------------------------------

      do j = jlo-1, jhi
      do i = ilo-1, ihi
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
  
         xcl = -0.5d0
         ycl =  0.d0

         xcr =  0.5d0
         ycr =  0.d0

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

         xdm = 0.5d0 * (xdr + xdl)
         ydm = 0.5d0 * (ydr + ydl)

         ! Intersection points

         xil = xcl
         yil = (xcl*(ydm-ydl) + xdm*ydl - xdl*ydm) / (xdm - xdl)
         
         xir = xcr
         yir = (xcr*(ydr-ydm) - xdm*ydr + xdr*ydm) / (xdr - xdm) 
         
         md = (ydr - ydl) / (xdr - xdl)
         
         if (abs(md) > puny) then
            xic = xdl - ydl/md
         else
            xic = 0.d0
         endif
         yic = 0.d0

         xicl = xic
         yicl = yic
         xicr = xic
         yicr = yic

    !-------------------------------------------------------------------
    ! Locate triangles in TL cell (NW for north edge, NE for east edge)
    ! and BL cell (W for north edge, N for east edge).
    !-------------------------------------------------------------------

         if (yil > 0.d0 .and. xdl < xcl .and. ydl >= 0.d0) then

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

         elseif (yil < 0.d0 .and. xdl < xcl .and. ydl < 0.d0) then

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

         elseif (yil < 0.d0 .and. xdl < xcl .and. ydl >= 0.d0) then

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

         elseif (yil > 0.d0 .and. xdl < xcl .and. ydl < 0.d0) then

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

         if (yir > 0.d0 .and. xdr >= xcr .and. ydr >= 0.d0) then

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

         elseif (yir < 0.d0 .and. xdr >= xcr .and. ydr < 0.d0) then

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

         elseif (yir < 0.d0 .and. xdr >= xcr  .and. ydr >= 0.d0) then 

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

         elseif (yir > 0.d0 .and. xdr >= xcr .and. ydr < 0.d0) then 

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
    ! For prescribed_area = T, shift the midpoint so that the departure
    ! region has the prescribed area
    !-------------------------------------------------------------------

         if (prescribed_area) then

            ! Sum the areas of the left and right triangles.
            ! Note that yp(i,j,1,ng) = 0 for all triangles, so we can
            !  drop those terms from the area formula.

            ng = 1
            area1 = 0.5d0 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
                               yp(i,j,3,ng)                   &
                            -  yp(i,j,2,ng) *                 &
                              (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
                            * areafact(i,j,ng) 

            ng = 2
            area2 = 0.5d0 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
                               yp(i,j,3,ng)                   &
                            -  yp(i,j,2,ng) *                 &
                              (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
                             * areafact(i,j,ng) 

            ng = 3
            area3 = 0.5d0 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
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

            if (ydl*ydr >= 0.d0) then   ! Both DPs lie on same side of x-axis

               ! compute required area of central departure region
               area_c  = edgearea(i,j) - area1 - area2 - area3

               ! shift midpoint so that the area of remaining triangles = area_c
               w1 = 2.d0*area_c/areafac_c(i,j)    &
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
                  xicl = 0.d0
               endif
               yicl = 0.d0

               if (abs(mdr) > puny) then
                  xicr = xdr - ydr/mdr
               else
                  xicr = 0.d0
               endif
               yicr = 0.d0

            elseif (xic < 0.d0) then  ! fix ICL = IC

               xicl = xic
               yicl = yic

               ! compute midpoint between ICL and DR
               xdm = 0.5d0 * (xdr + xicl)
               ydm = 0.5d0 *  ydr

               ! compute area of triangle adjacent to left corner 
               area4 = 0.5d0 * (xcl - xic) * ydl * areafac_l(i,j)
               area_c  = edgearea(i,j) - area1 - area2 - area3 - area4

               ! shift midpoint so that area of remaining triangles = area_c
               w1 = 2.d0*area_c/areafac_c(i,j) + (xcr-xic)*ydr
               w2 = (xdr-xic)**2 + ydr**2
               w1 = w1/w2
               xdm = xdm + ydr*w1
               ydm = ydm - (xdr - xic) * w1

               ! compute ICR
               mdr = (ydr - ydm) / (xdr - xdm)
               if (abs(mdr) > puny) then
                  xicr = xdr - ydr/mdr
               else
                  xicr = 0.d0
               endif
               yicr = 0.d0

            elseif (xic >= 0.d0) then  ! fix ICR = IR

               xicr = xic
               yicr = yic

               ! compute midpoint between ICR and DL 
               xdm = 0.5d0 * (xicr + xdl)
               ydm = 0.5d0 *  ydl

               area4 = 0.5d0 * (xic - xcr) * ydr * areafac_r(i,j)
               area_c  = edgearea(i,j) - area1 - area2 - area3 - area4

               ! shift midpoint so that area of remaining triangles = area_c
               w1 = 2.d0*area_c/areafac_c(i,j) + (xic-xcl)*ydl
               w2 = (xic-xdl)**2 + ydl**2
               w1 = w1/w2
               xdm = xdm - ydl*w1
               ydm = ydm - (xic - xdl) * w1

               ! compute ICL

               mdl = (ydm - ydl) / (xdm - xdl)
               if (abs(mdl) > puny) then
                  xicl = xdl - ydl/mdl
               else
                  xicl = 0.d0
               endif
               yicl = 0.d0

            endif   ! ydl*ydr >= 0.d0

         endif  ! prescribed_area

    !-------------------------------------------------------------------
    ! Locate triangles in BC cell (H for both north and east edges) 
    ! and TC cell (N for north edge and E for east edge).
    !-------------------------------------------------------------------

    ! Start with cases where both DPs lie in the same grid cell

         if (ydl >= 0.d0 .and. ydr >= 0.d0 .and. ydm >= 0.d0) then

         ! T1.d0a (group 4)

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

         elseif (ydl >= 0.d0 .and. ydr >= 0.d0 .and. ydm < 0.d0) then  ! rare

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

         elseif (ydl < 0.d0 .and. ydr < 0.d0 .and. ydm < 0.d0) then

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

         elseif (ydl < 0.d0 .and. ydr < 0.d0 .and. ydm >= 0.d0) then  ! rare

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

         elseif (ydl >= 0.d0 .and. ydr < 0.d0 .and. xic >= 0.d0  &
                                          .and. ydm >= 0.d0) then

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

         elseif (ydl >= 0.d0 .and. ydr < 0.d0 .and. xic >= 0.d0  &
                                          .and. ydm < 0.d0 ) then  ! less common

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

         elseif (ydl >= 0.d0 .and. ydr < 0.d0 .and. xic < 0.d0   &
                                          .and. ydm < 0.d0) then

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

         elseif (ydl >= 0.d0 .and. ydr < 0.d0 .and. xic <  0.d0  &
                                          .and. ydm >= 0.d0) then  ! less common

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

         elseif (ydl < 0.d0 .and. ydr >= 0.d0 .and. xic <  0.d0  &
                                          .and. ydm >= 0.d0) then

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

         elseif (ydl < 0.d0 .and. ydr >= 0.d0 .and. xic < 0.d0  &
                                          .and. ydm < 0.d0) then ! less common

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

         elseif (ydl < 0.d0 .and. ydr >= 0.d0 .and. xic >= 0.d0  &
                                          .and. ydm <  0.d0) then

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

         elseif (ydl < 0.d0 .and. ydr >= 0.d0 .and. xic >= 0.d0   &
                                          .and. ydm >= 0.d0) then  ! less common

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
    ! If prescribed_area = T, this area is compared to edgearea as a bug check.
    ! If prescribed_area = F, this area is passed as an output array.
    !-------------------------------------------------------------------

      areasum(:,:) = 0.d0

      do ng = 1, ngroups
         icells(ng) = 0

         do ij = 1, icellsd
            i = indxid(ij)
            j = indxjd(ij)

            triarea(i,j,ng) = 0.5d0 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
                                     (yp(i,j,3,ng)-yp(i,j,1,ng))   &
                                   - (yp(i,j,2,ng)-yp(i,j,1,ng)) *   &
                                     (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
                                   * areafact(i,j,ng) 

            if (abs(triarea(i,j,ng)) < 1.e-16*areafac_c(i,j)) then
               triarea(i,j,ng) = 0.d0
            else
               icells(ng) = icells(ng) + 1 
               ic = icells(ng)
               indxi(ic,ng) = i
               indxj(ic,ng) = j
            endif

            areasum(i,j) = areasum(i,j) + triarea(i,j,ng)

         enddo                  ! ij
      enddo                     ! ng

      if (prescribed_area) then
       if (bugcheck) then   ! set bugcheck = F to speed up code
         do ij = 1, icellsd
            i = indxid(ij)
            j = indxjd(ij)
            if (abs(areasum(i,j) - edgearea(i,j)) > 1.e-13*areafac_c(i,j)) then
               print*, ''
               print*, 'Areas do not add up: i, j, edge =',   &
                        i, j, trim(edge)
               print*, 'edgearea =', edgearea(i,j)
               print*, 'areasum =', areasum(i,j)
               print*, 'areafac_c =', areafac_c(i,j)
               print*, ''
               print*, 'Triangle areas:'
               do ng = 1, ngroups   
                  if (abs(triarea(i,j,ng)) > 1.e-16*abs(areafact(i,j,ng))) then
                     print*, ng, triarea(i,j,ng)
                  endif
               enddo
            endif
         enddo
       endif          ! bugcheck

      else            ! prescribed_area = F
         do ij = 1, icellsd
            i = indxid(ij)
            j = indxjd(ij)
            edgearea(i,j) = areasum(i,j)
         enddo
      endif     ! prescribed_area

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
                  xp(i,j,nv,ng) = xp(i,j,nv,ng) - 1.d0*ishift
                  yp(i,j,nv,ng) = yp(i,j,nv,ng) + 0.5d0 - 1.d0*jshift
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
                  xp(i,j,nv,ng) =  yp(i,j,nv,ng) + 0.5d0 - 1.d0*ishift
                  yp(i,j,nv,ng) = -w1 - 1.d0*jshift
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
                  if (abs(xp(i,j,nv,ng)) > 0.5d0+puny) then
                     print*, ''
                     print*, 'WARNING: xp =', xp(i,j,nv,ng)
                     print*, 'i, j, ng, nv =', i, j, ng, nv
!                     print*, 'yil,xdl,xcl,ydl=',yil,xdl,xcl,ydl
!                     print*, 'yir,xdr,xcr,ydr=',yir,xdr,xcr,ydr
!                     print*, 'ydm=',ydm
!                      stop
                  endif
                  if (abs(yp(i,j,nv,ng)) > 0.5d0+puny) then
                     print*, ''
                     print*, 'WARNING: yp =', yp(i,j,nv,ng)
                     print*, 'i, j, ng, nv =', i, j, ng, nv
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
                                       icells,                    &
                                       indxi,          indxj,     &
                                       xp,             yp,        &
                                       integral_order)
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
      integer, intent(in) ::   &
           nx_block, ny_block  ! block dimensions

      integer, dimension (ngroups), intent(in) ::     &
           icells              ! number of cells where triarea > puny

      integer, dimension (nx_block*ny_block,ngroups),     &
           intent(in) ::     &
           indxi ,&! compressed index in i-direction
           indxj   ! compressed index in j-direction

      real(dp), intent(inout),   &
           dimension (nx_block, ny_block, 0:nvert, ngroups) ::   &
           xp, yp          ! coordinates of triangle points

      integer, intent(in) ::   &
           integral_order  ! 1 = linear, 2 = quadratic
!
!EOP
!
      integer ::   &
           i, j, ij          ,&! horizontal indices
           ng                  ! triangle index

      if (integral_order == 1) then ! linear (1-point formula)

         do ng = 1, ngroups
         do ij = 1, icells(ng)
            i = indxi(ij,ng)
            j = indxj(ij,ng)

            ! coordinates of midpoint
            xp(i,j,0,ng) = (xp(i,j,1,ng) + xp(i,j,2,ng) + xp(i,j,3,ng)) / 3.d0
            yp(i,j,0,ng) = (yp(i,j,1,ng) + yp(i,j,2,ng) + yp(i,j,3,ng)) / 3.d0

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
            xp(i,j,0,ng) = (xp(i,j,1,ng) + xp(i,j,2,ng) + xp(i,j,3,ng)) / 3.d0
            yp(i,j,0,ng) = (yp(i,j,1,ng) + yp(i,j,2,ng) + yp(i,j,3,ng)) / 3.d0

            ! coordinates of the 3 points needed for integrals

            xp(i,j,1,ng) = 0.5d0*xp(i,j,1,ng) + 0.5d0*xp(i,j,0,ng)
            yp(i,j,1,ng) = 0.5d0*yp(i,j,1,ng) + 0.5d0*yp(i,j,0,ng)

            xp(i,j,2,ng) = 0.5d0*xp(i,j,2,ng) + 0.5d0*xp(i,j,0,ng)
            yp(i,j,2,ng) = 0.5d0*yp(i,j,2,ng) + 0.5d0*yp(i,j,0,ng)

            xp(i,j,3,ng) = 0.5d0*xp(i,j,3,ng) + 0.5d0*xp(i,j,0,ng)
            yp(i,j,3,ng) = 0.5d0*yp(i,j,3,ng) + 0.5d0*yp(i,j,0,ng)

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
      subroutine transport_integrals (nx_block,       ny_block,       &
                                      ntracer,        icells,         &
                                      indxi,          indxj,          &
                                      triarea,        integral_order, &
                                      iflux,          jflux,          &
                                      xp,             yp,             &
                                      mc,             mx,             &
                                      my,             mflx,           &
                                      tc,             tx,             &
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
      integer, intent(in) ::   &
           nx_block, ny_block  ,&! block dimensions
           ntracer               ! number of tracers in use

      integer, dimension (ngroups), intent(in) ::     &
           icells           ! number of cells where triarea > puny

      integer, dimension (nx_block*ny_block,ngroups),     &
           intent(in) ::     &
           indxi ,&! compressed index in i-direction
           indxj   ! compressed index in j-direction

      real(dp), intent(in),   &
           dimension (nx_block, ny_block, 0:nvert, ngroups) ::   &
           xp, yp           ! coordinates of triangle points

      real(dp), intent(in),   &
           dimension (nx_block, ny_block, ngroups) ::   &
           triarea          ! triangle area

      integer, intent(in) ::   &
           integral_order  ! 1 = linear, 2 = quadratic

      integer, intent(in),   &
           dimension (nx_block, ny_block, ngroups) ::   &
           iflux     ,&
           jflux

      real(dp), intent(in),   &
           dimension (nx_block, ny_block) ::   &
           mc, mx, my

      real(dp), intent(out),   &
           dimension (nx_block, ny_block) ::   &
           mflx

      real(dp), intent(in),   &
           dimension (nx_block, ny_block, ntracer), optional ::   &
           tc, tx, ty

      real(dp), intent(out),   &
           dimension (nx_block, ny_block, ntracer), optional ::   &
           mtflx
!
!EOP
!
      integer ::   &
           i, j, ij      ,&! horizontal indices of edge
           i2, j2        ,&! horizontal indices of cell contributing transport
           ng            ,&! triangle index
           nt, nt1       ,&! tracer indices
           ilo,ihi,jlo,jhi ! beginning and end of physical domain

      real(dp) ::   &
           m0, m1, m2, m3         ,&! mass field at internal points
           w0, w1, w2, w3           ! work variables

      real(dp), dimension (nx_block, ny_block) ::   &
           msum, mxsum, mysum     ,&! sum of mass, mass*x, and mass*y
           mxxsum, mxysum, myysum   ! sum of mass*x*x, mass*x*y, mass*y*y

      real(dp), dimension (nx_block, ny_block, ntracer) ::   &
           mtsum                    ! sum of mass*tracer

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

      mflx(:,:) = 0.d0
      if (present(mtflx)) then
         do nt = 1, ntracer
           mtflx(:,:,nt) = 0.d0
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
               m1 = (mc(i2,j2) + xp(i,j,1,ng)*mx(i2,j2)   &
                               + yp(i,j,1,ng)*my(i2,j2)) / 3.d0
               m2 = (mc(i2,j2) + xp(i,j,2,ng)*mx(i2,j2)   &
                               + yp(i,j,2,ng)*my(i2,j2)) / 3.d0
               m3 = (mc(i2,j2) + xp(i,j,3,ng)*mx(i2,j2)   &
                               + yp(i,j,3,ng)*my(i2,j2)) / 3.d0
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

         endif                  ! integral_order

         ! mass * tracer transports

         if (present(mtflx)) then

            do nt = 1, ntracer

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

               enddo         ! ij

            enddo               ! ntracer
        endif                   ! present(mtflx)
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
                                ntracer,                 &
                                tarear,      l_stop,     &
                                mflxe,       mflxn,      &
                                mass,                    &
                                mtflxe,      mtflxn,     &
                                trcr)
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
      integer, intent(in) ::   &
         nx_block, ny_block,&! block dimensions
         ilo,ihi,jlo,jhi   ,&! beginning and end of physical domain
         ntracer             ! number of tracers in use

      real(dp), dimension (nx_block, ny_block), intent(in) ::   &
         mflxe, mflxn     ! mass transport across east and north cell edges

      real(dp), intent(in) ::   &
         tarear           ! 1/tarea

      real(dp), dimension (nx_block, ny_block),   &
         intent(inout) ::   &
         mass             ! mass field (mean)

      real(dp), dimension (nx_block, ny_block, ntracer),  &
         intent(in), optional ::   &
         mtflxe, mtflxn   ! mass*tracer transport across E and N cell edges

      real(dp), dimension (nx_block, ny_block, ntracer),  &
         intent(inout), optional ::   &
         trcr             ! tracer fields

      logical, intent(inout) ::   &
         l_stop           ! if true, abort on return
!
!EOP
!
      integer ::   &
         i, j           ,&! horizontal indices
         nt, nt1, nt2     ! tracer indices

      real(dp), dimension(nx_block,ny_block,ntracer) ::   &
         mtold            ! old mass*tracer

      real(dp) ::   &
         w1, w2           ! work variables

      integer, dimension(nx_block*ny_block) ::   &
         indxi          ,&! compressed indices in i and j directions
         indxj

      integer ::   &
         icells         ,&! number of cells with mass > 0.
         ij               ! combined i/j horizontal index

      character(len=100) :: message

      integer ::   &
         istop, jstop     ! indices of grid cell where model aborts 

    !-------------------------------------------------------------------
    ! Save starting values of mass*tracer
    !-------------------------------------------------------------------

      if (present(trcr)) then
         do nt = 1, ntracer
           do j = jlo, jhi
            do i = ilo, ihi
               mtold(i,j,nt) = mass(i,j) * trcr(i,j,nt)
            enddo            ! i
            enddo              ! j
         enddo                  ! nt
      endif                     ! present(trcr)

    !-------------------------------------------------------------------
    ! Update mass field
    !-------------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi

         w1 = mflxe(i,j) - mflxe(i-1,j)   &
            + mflxn(i,j) - mflxn(i,j-1)
         mass(i,j) = mass(i,j) - w1*tarear

         if (mass(i,j) < -puny) then    ! abort with negative value
            l_stop = .true.
            istop = i
            jstop = j
         elseif (mass(i,j) < 0.d0) then   ! set to zero
            mass(i,j) = 0.d0
         endif

      enddo
      enddo

!WHL - Test the diagnostics
!TODO - Write error message cleanly to log file.
!       For now, just print out an error message.

      if (l_stop) then
         i = istop
         j = jstop
         w1 = mflxe(i,j) - mflxe(i-1,j)   &
            + mflxn(i,j) - mflxn(i,j-1)
!         write (message,*) 'Process:',this_rank
!         call write_log(message)
!         write (message,*) 'Remap, negative ice thickness, i, j =', i, j
!         call write_log(message)    
!         write (message,*) 'Old thickness =', mass(i,j) + w1*tarear
!         call write_log(message)    
!         write (message,*) 'New thickness =', mass(i,j)
!         call write_log(message)    
!         write (message,*) 'Net transport =', -w1*tarear
!         call write_log(message)    
         write (6,*) 'Process:',this_rank
         write (6,*) 'Remap, negative ice thickness, i, j =', i, j
         write (6,*) 'Old thickness =', mass(i,j) + w1*tarear
         write (6,*) 'New thickness =', mass(i,j)
         write (6,*) 'Net transport =', -w1*tarear
         return
      endif

    !-------------------------------------------------------------------
    ! Update tracers
    !-------------------------------------------------------------------

      if (present(trcr)) then

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (mass(i,j) > 0.d0) then ! grid cells with positive areas
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo                  ! i
         enddo                  ! j

         do nt = 1, ntracer

           do j = jlo, jhi
            do i = ilo, ihi
               trcr(i,j,nt) = 0.d0
            enddo
            enddo

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               w1  = mtflxe(i,j,nt) - mtflxe(i-1,j,nt)   &
                   + mtflxn(i,j,nt) - mtflxn(i,j-1,nt)
               trcr(i,j,nt) = (mtold(i,j,nt) - w1*tarear)   &
                             / mass(i,j)
            enddo            ! ij

         enddo                  ! nt
      endif                     ! present(trcr)

      end subroutine update_fields

!=======================================================================

      end module glissade_remap

!=======================================================================
