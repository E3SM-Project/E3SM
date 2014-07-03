!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_transport.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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
! !MODULE: glissade_transport - driver for transport of mass and tracers
!
! !DESCRIPTION:
!
! Drivers for remapping and upwind ice transport
!
! author: William H. Lipscomb, LANL 
!
! This version was created from ice_transport_driver in CICE,
!  revision 313, 6 Jan. 2011.
! The repository is here: http://oceans11.lanl.gov/svn/CICE

! !INTERFACE:

      module glissade_transport
!
! !USES:
      use glimmer_global, only: dp
      use glimmer_log
      use glissade_remap, only: glissade_horizontal_remap, make_remap_mask, puny

      use parallel 

!
!EOP
!
      implicit none
      save
      private
      public :: glissade_transport_driver

      logical, parameter ::  &
         prescribed_area = .false.  ! if true, prescribe the area fluxed across each edge

!TODO - Code uses Protex documenting.  Revise for doxygen

!WHL = debug
!!      integer, parameter :: idiag = 31, jdiag = 31


!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: glissade_transport_driver - driver for mass and tracer transport
!
! !INTERFACE:
!
      subroutine glissade_transport_driver(dt,                   &
                                           dx,       dy,         &
                                           nx,       ny,         &
                                           nlyr,     sigma,      &
                                           nhalo,    ntracer,    &
                                           uvel,     vvel,       &
                                           thck,                 &
                                           acab,     bmlt,       &
                                           temp,     age,        &
                                           upwind_transport_in)

!
! !DESCRIPTION:
!
! This subroutine solves the transport equations for one timestep
! using the conservative remapping scheme developed by John Dukowicz
! and John Baumgardner and modified for sea ice by William
! Lipscomb and Elizabeth Hunke.
!
! This scheme preserves monotonicity of ice area and tracers.  That is,
! it does not produce new extrema.  It is second-order accurate in space,
! except where gradients are limited to preserve monotonicity. 
!
! Optionally, the remapping scheme can be replaced with a simple
! first-order upwind scheme.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!

      real(dp), intent(in) ::     &
         dt,                   &! time step (s)
         dx, dy                 ! gridcell dimensions (m)
                                ! (cells assumed to be rectangular)

      !TODO - Use nhalo in parallel module instead of passing in
      !     - Declare ntracer somewhere instead of passing in?

      integer, intent(in) ::   &
         nx, ny,               &! horizontal array size
         nlyr,                 &! number of vertical layers
         ntracer,              &! number of tracers (should be >= 1)
         nhalo                  ! number of rows of halo cells

      real(dp), intent(in), dimension(nlyr+1) :: &
         sigma                  ! layer interfaces in sigma coordinates
                                ! top sfc = 0, bottom sfc = 1

      !Note dimensions of uvel and vvel
      real(dp), intent(in), dimension(nlyr+1,nx-1,ny-1) :: &
         uvel, vvel             ! horizontal velocity components (m/s)
                                ! (defined at horiz cell corners, vertical interfaces)

      real(dp), intent(inout), dimension(nx,ny) :: &
         thck                   ! ice thickness (m)
                                ! (defined at horiz cell centers)

      real(dp), intent(in), dimension(nx,ny) :: &
         acab                   ! surface mass balance (m/s)
                                ! (defined at horiz cell centers)

      real(dp), intent(in), dimension(nx,ny) :: &
         bmlt                   ! basal melt rate (m/s)
                                ! positive for melting, negative for freeze-on
                                ! (defined at horiz cell centers)

      ! Note vertical dimension of temperature
      real(dp), intent(inout), dimension(0:nlyr+1,nx,ny), optional :: &
         temp                   ! ice temperature
                                ! (defined at horiz cell centers, vertical layer midpts)

      real(dp), intent(inout), dimension(nlyr,nx,ny), optional :: &
         age                    ! ice age

      logical, intent(in), optional ::  &
         upwind_transport_in    ! if true, do first-order upwind transport
! 
!EOP
!
      ! local variables

      integer ::     &
         i, j, k         ,&! cell indices
         ilo,ihi,jlo,jhi ,&! beginning and end of physical domain
         nt                ! tracer index

      real(dp), dimension (nx,ny) ::     &
         thck_mask         ! = 1. if ice is present, = 0. otherwise

      real(dp), dimension (nx-1,ny-1) ::     &
         uvel_layer      ,&! uvel averaged to layer midpoint (m/s)
         vvel_layer        ! vvel averaged to layer midpoint (m/s)

      real(dp), dimension (nx,ny,nlyr) ::     &
         thck_layer        ! ice layer thickness (m)

      real(dp), dimension (nx,ny,ntracer,nlyr) ::     &
         tracer            ! tracer values

      real(dp), dimension (nx,ny) ::     &
         tsfc,            &! surface temperature, temp(0,:,:)
         tbed              ! basal temperature, temp(nlyr+1,:,:)

      integer ::     &
         icells            ! number of cells with ice

      integer, dimension(nx*ny) ::     &
         indxi, indxj      ! compressed i/j indices

      real(dp) ::   &
         new_thck,        &! temporary ice thickness
         new_bmlt          ! bmlt adjusted for residual ablation

      real(dp), dimension (nlyr) ::     &
         new_thck_layer    ! temporary ice layer thickness

    !-------------------------------------------------------------------
    ! If prescribed_area is true, the area of each departure region is
    !  computed in advance (e.g., by taking the divergence of the 
    !  velocity field and passed to locate_triangles.  The departure 
    !  regions are adjusted to obtain the desired area.
    ! If false, edgearea is computed in locate_triangles and passed out.
    !-------------------------------------------------------------------

      real(dp), dimension(nx,ny) ::   &
         edgearea_e     ,&! area of departure regions for east edges
         edgearea_n       ! area of departure regions for north edges

      logical, parameter ::     &
         conservation_check = .true. ! if true, check global conservation

      !TODO - separate arrays for mass and tracers?
      real(dp), dimension(0:ntracer) ::     &
         mtsum_init,     &! initial global ice mass and global ice mass*tracer
         mtsum_final      ! final global ice mass and global ice mass*tracer

      real(dp) ::    &
         sfc_accum, sfc_ablat,  &! accumulation and ablation at upper surface (m)
         bed_accum, bed_ablat    ! freeze-on and melting at lower surface (m)

      logical ::     &
         l_stop            ! if true, abort the model

      character(len=100) :: message

      real(dp), dimension (:,:,:), allocatable :: &
         worku            ! work array

      real(dp), dimension(nx,ny) ::      &
         uee, vnn            ! cell edge velocities for upwind transport

      logical ::     &
         upwind_transport    ! if true, do first-order upwind transport

      if (present(upwind_transport_in)) then
         upwind_transport = upwind_transport_in
      else
         upwind_transport = .false.
      endif

!---!-------------------------------------------------------------------
!---! Prepare for remapping.
!---! Initialize, update halo cells, fill tracer arrays.
!---!-------------------------------------------------------------------

!Note: (ilo,ihi) and (jlo,jhi) are the lower and upper bounds of the local domain
! (i.e., grid cells owned by this processor).

!TODO - Replace ilo with 1+nhalo, ihi with nx-nhalo, etc.?

      ilo = 1 + nhalo
      ihi = nx - nhalo
      jlo = 1 + nhalo
      jhi = ny - nhalo

      l_stop = .false.

    !-------------------------------------------------------------------
    ! Fill thickness and tracer arrays.
    ! Assume that temperature (if present) is tracer 1, and age (if present)
    !  is tracer 2.  Add more tracers as needed.
    ! If no tracers are present, then only the ice thickness is transported.
    !  In this case we define a dummy tracer array, since glissade_horizontal_remap
    !  requires that a tracer array be passed in.
    !-------------------------------------------------------------------

      do k = 1, nlyr

         thck_layer(:,:,k) = thck(:,:) * (sigma(k+1) - sigma(k))

         if (present(temp)) then
            tracer(:,:,1,k) = temp(k,:,:)
         else
            tracer(:,:,1,k) = 0.d0    ! dummy array
         endif

         if (present(age) .and. ntracer >= 2) tracer(:,:,2,k) = age(k,:,:)

         !TODO - Other tracer fields could be added here

      enddo

    !-------------------------------------------------------------------
    ! Set surface and basal temperatures.
    !-------------------------------------------------------------------

      if (present(temp)) then
         tsfc(:,:) = temp(0,:,:)
         tbed(:,:) = temp(nlyr+1,:,:)
      else
         tsfc(:,:) = 0.d0
         tbed(:,:) = 0.d0
      endif

    !TODO - Put this conservation check code in a subroutine?
    !-------------------------------------------------------------------
    ! Compute initial values of globally conserved quantities (optional)
    !-------------------------------------------------------------------

      if (conservation_check) then

         ! Compute initial values of globally conserved quantities,
         !  accounting for surface accumulation.
         ! Assume gridcells of equal area, ice of uniform density.

         mtsum_init(:) = 0.0_dp

         do j = jlo, jhi
         do i = ilo, ihi

            !------------------------------------------------------------------
            ! accumulate ice thickness (proportional to mass) and thickness*tracers
            !------------------------------------------------------------------

            mtsum_init(0) = mtsum_init(0) + thck(i,j)   
            do nt = 1, ntracer
               mtsum_init(nt) =  mtsum_init(nt) + sum(tracer(i,j,nt,:) * thck_layer(i,j,:))
            enddo

            !------------------------------------------------------------------
            ! account for accumulation/ablation at top surface
            ! assume acab > 0 for accumulation, < 0 for ablation
            !------------------------------------------------------------------

            new_thck_layer(:) = thck_layer(i,j,:)   ! keep track of new layer thicknesses as we handle acab,
                                                    ! so we can deal correctly with bmlt below

            if (acab(i,j) > 0.d0) then  ! net accumulation

               sfc_accum = acab(i,j)*dt
               new_thck_layer(1) = new_thck_layer(1) + sfc_accum

               ! add mass
               mtsum_init(0) = mtsum_init(0) + sfc_accum

               ! add mass*temperature
               mtsum_init(1) = mtsum_init(1) + tsfc(i,j) * sfc_accum

               ! assume age of new ice = 0, so do not need to increment mtsum_init(2)

               !TODO - Add other tracers here, as needed

            elseif (acab(i,j) < 0.d0) then  ! net ablation
             
               sfc_ablat = -acab(i,j)*dt

               if (thck(i,j) > sfc_ablat) then

                  ! subtract mass
                  mtsum_init(0) = mtsum_init(0) - sfc_ablat

                  ! subtract mass*tracers
                  do k = 1, nlyr      
                     if (sfc_ablat > thck_layer(i,j,k)) then  ! remove entire layer
                        do nt = 1, ntracer
                           mtsum_init(nt) = mtsum_init(nt) - tracer(i,j,nt,k)*thck_layer(i,j,k)
                        enddo
                        new_thck_layer(k) = 0.d0
                        sfc_ablat = sfc_ablat - thck_layer(i,j,k)
                     else                                     ! remove part of layer
                        do nt = 1, ntracer
                           mtsum_init(nt) = mtsum_init(nt) - tracer(i,j,nt,k)*sfc_ablat
                        enddo
                        new_thck_layer(k) = new_thck_layer(k) - sfc_ablat
                        sfc_ablat = 0.d0
                        exit
                     endif
                  enddo

               else    ! thck < sfc_ablat, so entire column will be removed

                  mtsum_init(0) = mtsum_init(0) - thck(i,j)
                  do k = 1, nlyr
                     do nt = 1, ntracer
                        mtsum_init(nt) = mtsum_init(nt) - tracer(i,j,nt,k)*thck_layer(i,j,k)
                     enddo
                     new_thck_layer(k) = 0.d0
                  enddo

                  ! There is some residual melt potential which we ignore for now.
                  ! May have to revisit this for coupled runs. 
                  sfc_ablat = sfc_ablat - sum(thck_layer(i,j,:))
                                                                   
               endif   ! thck > sfc_ablat
 
            endif  ! acab > 0

            !------------------------------------------------------------------
            ! account for melting/freeze-on at bed
            ! assume bmlt is > 0 for melting, < 0 for freezing
            !------------------------------------------------------------------

            ! new ice thickness adjusted for acab
            new_thck = sum(new_thck_layer(:))
 
            if (bmlt(i,j) < 0.d0) then  ! freeze-on

               bed_accum = -bmlt(i,j)*dt

               ! add mass
               mtsum_init(0) = mtsum_init(0) + bed_accum

               ! add mass*temperature
               mtsum_init(1) = mtsum_init(1) + tbed(i,j) * bed_accum

               ! assume age of new ice = 0, so do not need to increment mtsum_init(2)

               !TODO - Add other tracers here, as needed

            elseif (bmlt(i,j) > 0.d0) then  ! melting
             
               bed_ablat = bmlt(i,j)*dt

               if (new_thck > bed_ablat) then

                  ! subtract mass
                  mtsum_init(0) = mtsum_init(0) - bed_ablat

                  ! subtract mass*tracers
                  do k = nlyr, 1, -1      
                     if (bed_ablat > new_thck_layer(k)) then  ! remove entire layer
                        do nt = 1, ntracer
                           mtsum_init(nt) = mtsum_init(nt) - tracer(i,j,nt,k)*new_thck_layer(k)
                        enddo
                        bed_ablat = bed_ablat - new_thck_layer(k)
                     else                                     ! remove part of layer
                        do nt = 1, ntracer
                           mtsum_init(nt) = mtsum_init(nt) - tracer(i,j,nt,k)*bed_ablat
                        enddo
                        bed_ablat = 0.d0
                        exit
                     endif
                  enddo

               else    ! new_thck < bed_ablat, so entire column will be removed

                  mtsum_init(0) = mtsum_init(0) - new_thck
                  do k = 1, nlyr
                     do nt = 1, ntracer
                        mtsum_init(nt) = mtsum_init(nt) - tracer(i,j,nt,k)*new_thck_layer(k)
                     enddo
                  enddo

                  ! There is some residual melt potential which we ignore for now.
                  ! May have to revisit this for coupled runs.
                  bed_ablat = bed_ablat - new_thck

               endif   ! thck > bed_ablat
 
            endif  ! bmlt < 0

         enddo       ! i
         enddo       ! j

         call global_sum(mtsum_init)

      endif                     ! conservation_check

      !TODO - Test upwind transport option

      if (upwind_transport) then

         allocate (worku(nx,ny,0:ntracer))

         do k = 1, nlyr

            ! Average corner velocities at layer interfaces to 
            ! edge velocities at layer midpoints.
      
            do j = jlo, jhi
            do i = ilo-1, ihi   ! include west edge of local domain
               uee(i,j) = 0.5d0 * (uvel(k,  i,j) + uvel(k,  i,j-1)    &
                                 + uvel(k+1,i,j) + uvel(k+1,i,j-1))
            enddo
            enddo


            do j = jlo-1, jhi   ! include south edge of local domain
            do i = ilo, ihi
               vnn(i,j) = 0.5d0 * (vvel(k,  i,j) + vvel(k,  i-1,j)    &
                                 + vvel(k+1,i,j) + vvel(k+1,i-1,j))
            enddo
            enddo

            ! Fill work array for transport

            worku(:,:,0) = thck_layer(:,:,k)
            do nt = 1, ntracer
               worku(:,:,nt) = thck_layer(:,:,k) * tracer(:,:,nt,k)
            enddo

            !-----------------------------------------------------------------
            ! Upwind transport
            !-----------------------------------------------------------------

            do nt = 0, ntracer
               call upwind_field (nx,             ny,                  &
                                  ilo, ihi,       jlo, jhi,            &
                                  dx,             dy,                  &
                                  dt,             worku(:,:,nt),       &
                                  uee(:,:),       vnn    (:,:))

            enddo   ! ntracer

            ! Recompute thickness and tracers

            thck_layer(:,:,k) = worku(:,:,0)
            do nt = 1, ntracer
               do j = jlo, jhi
               do i = ilo, ihi
                  if (thck_layer(i,j,k) > puny) then
                     tracer(i,j,nt,k) = worku(i,j,nt) / thck_layer(i,j,k)
                  else
                     tracer(i,j,nt,k) = 0.d0
                  endif
               enddo   ! i
               enddo   ! j
            enddo      ! ntracer

         enddo         ! nlyr

         deallocate (worku)

      else    ! remapping transport

      !-------------------------------------------------------------------
      ! Define a mask: = 1 where ice is present, = 0 otherwise         
      ! The mask is used to prevent tracer values in cells without ice from
      !  being used to compute tracer gradients.
      !-------------------------------------------------------------------

         call make_remap_mask (nx,           ny,                 &
                               ilo, ihi,     jlo, jhi,           &
                               nhalo,        icells,             &
                               indxi(:),     indxj(:),           &
                               thck(:,:),    thck_mask(:,:))

      !-------------------------------------------------------------------
      ! Remap ice thickness and tracers; loop over layers
      !-------------------------------------------------------------------

         do k = 1, nlyr

            ! Average velocities to the midpoint of the layer

            uvel_layer(:,:) = 0.5d0 * (uvel(k,:,:) + uvel(k+1,:,:))
            vvel_layer(:,:) = 0.5d0 * (vvel(k,:,:) + vvel(k+1,:,:))

            edgearea_e(:,:) = 0.d0
            edgearea_n(:,:) = 0.d0

            ! If prescribed_area is true, compute edgearea by taking the divergence
            !  of the velocity field.

            if (prescribed_area) then

               do j = jlo, jhi
               do i = ilo-1, ihi
                  edgearea_e(i,j) = (uvel_layer(i,j) + uvel_layer(i,j-1))     &
                                    * 0.5d0 * dy * dt
               enddo
               enddo
 
               do j = jlo-1, jhi
               do i = ilo, ihi
                  edgearea_n(i,j) = (vvel_layer(i,j) + vvel_layer(i-1,j))     &
                                   * 0.5d0 * dx * dt
               enddo
               enddo

            endif

         !-------------------------------------------------------------------
         ! Main remapping routine: Step ice thickness and tracers forward in time.
         !-------------------------------------------------------------------

            call glissade_horizontal_remap (dt,                                  &
                                            dx,                dy,               &
                                            nx,                ny,               &
                                            ntracer,           nhalo,            &
                                            thck_mask(:,:),    icells,           &
                                            indxi(:),          indxj(:),         &
                                            uvel_layer(:,:),   vvel_layer(:,:),  &
                                            thck_layer(:,:,k), tracer(:,:,:,k),  &
                                            edgearea_e(:,:),   edgearea_n(:,:))

         enddo    ! nlyr

      endif       ! remapping v. upwind transport

      !-------------------------------------------------------------------
      ! Add the mass balance at the surface and bed.
      ! The reason to do this now rather than at the beginning of the
      !  subroutine is that the velocity was computed for the old geometry,
      !  before the addition or loss of new mass at the surface and bed.  
      ! TODO: Rethink this ordering if we move to implicit or semi-implicit
      !       timestepping, where the velocity depends on the new geometry.
      !
      ! We assume here that new ice arrives at the surface with the same 
      !  temperature as the surface.
      ! TODO: Make sure that this assumption is consistent with energy
      !       conservation for coupled simulations.
      !-------------------------------------------------------------------

      call glissade_add_smb(nx,       ny,         &
                            nlyr,     ntracer,    &
                            dt,                   &
                            thck_layer(:,:,:),    &
                            tracer(:,:,:,:),      &
                            acab(:,:),            &
                            tsfc(:,:),            &
                            bmlt(:,:),            &
                            tbed(:,:) )
   
      !-------------------------------------------------------------------
      ! Interpolate tracers back to sigma coordinates 
      !-------------------------------------------------------------------

      call glissade_vertical_remap(nx,                ny,      &
                                   nlyr,              ntracer, &
                                   sigma(:),                   &
                                   thck_layer(:,:,:),          &
                                   tracer(:,:,:,:) )

      !-------------------------------------------------------------------
      ! Halo updates for thickness and tracer arrays
      !
      ! Note: Cannot pass the full 3D array to parallel_halo, because this
      !       subroutine assumes that k is the first rather than third index. 
      !-------------------------------------------------------------------

      do k = 1, nlyr

         call parallel_halo(thck_layer(:,:,k))

         do nt = 1, ntracer
            call parallel_halo(tracer(:,:,nt,k))
         enddo

      enddo

      !-------------------------------------------------------------------
      ! Recompute thickness, temperature and other tracers.
      !-------------------------------------------------------------------

      thck(:,:) = 0.d0

      do k = 1, nlyr

         thck(:,:) = thck(:,:) + thck_layer(:,:,k)

         if (present(temp)) temp(k,:,:) = tracer(:,:,1,k)
         if (present(age) .and. ntracer >= 2) age(k,:,:) = tracer(:,:,2,k)
         !WHL - Could add more tracer fields here

      enddo

    !-------------------------------------------------------------------
    ! Compute final values of globally conserved quantities.
    ! Check global conservation of mass and mass*tracers.  (Optional)
    !-------------------------------------------------------------------

      if (conservation_check) then

         ! Compute final values of globally conserved quantities.
         ! Assume gridcells of equal area, ice of uniform density.

         mtsum_final(:) = 0.0_dp
         do j = jlo, jhi
         do i = ilo, ihi
            mtsum_final(0) = mtsum_final(0) + thck(i,j)  ! local for now
            do nt = 1, ntracer
               mtsum_final(nt) =  mtsum_final(nt) + sum(tracer(i,j,nt,:) * thck_layer(i,j,:))
            enddo
         enddo
         enddo
         call global_sum(mtsum_final)

         if (main_task) then
            call global_conservation (mtsum_init(0), mtsum_final(0),  &
                                      l_stop,        ntracer,         &
                                      mtsum_init(1:ntracer), mtsum_final(1:ntracer))
            if (l_stop) then
               write(message,*) 'CONSERVATION ERROR in glissade_transport; Aborting'
               call write_log(message,GM_FATAL)
!               call write_log(message,GM_DIAGNOSTIC)  ! uncomment for debugging
            endif
         endif

      endif                     ! conservation_check

      end subroutine glissade_transport_driver

!=======================================================================
!
!BOP
!
! !IROUTINE: global_conservation - check for changes in conserved quantities
!
! !INTERFACE:
!
      subroutine global_conservation (msum_init,  msum_final,     &
                                      l_stop,     ntracer,        &
                                      mtsum_init, mtsum_final)
!
! !DESCRIPTION:
!
! Check whether values of conserved quantities have changed.
! An error probably means that ghost cells are treated incorrectly.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real(dp), intent(in) ::     &
         msum_init   ,&! initial global ice mass 
         msum_final    ! final global ice mass

      logical, intent(inout) ::     &
         l_stop    ! if true, abort on return

      integer, intent(in) ::  &
         ntracer       ! number of tracers

      real(dp), dimension(ntracer), intent(in), optional :: &
         mtsum_init  ,&! initial global ice mass*tracer
         mtsum_final   ! final global ice mass*tracer

      character(len=100) :: message
!
!EOP
!
      integer ::     &
           nt            ! tracer index

      real(dp) ::     &
           diff          ! difference between initial and final values

      if (msum_init > puny) then
         diff = msum_final - msum_init
         if (abs(diff/msum_init) > puny) then
            l_stop = .true.
            write (message,*) 'glissade_transport: ice mass conservation error'
            call write_log(message)
            write (message,*) 'Initial global mass =', msum_init
            call write_log(message)
            write (message,*) 'Final global mass =', msum_final
            call write_log(message)
            write (message,*) 'Fractional error =', abs(diff)/msum_init
            call write_log(message)
         endif
      endif

      if (present(mtsum_init)) then
       do nt = 1, ntracer
         if (abs(mtsum_init(nt)) > puny) then
            diff = mtsum_final(nt) - mtsum_init(nt)
            if (abs(diff/mtsum_init(nt)) > puny) then
               l_stop = .true.
               write (message,*) 'glissade_transport: mass*tracer conservation error'
               call write_log(message)
               write (message,*) 'tracer index =', nt
               call write_log(message)
               write (message,*) 'Initial global mass*tracer =', mtsum_init(nt)
               call write_log(message)
               write (message,*) 'Final global mass*tracer =', mtsum_final(nt)
               call write_log(message)
               write (message,*) 'Fractional error =', abs(diff)/mtsum_init(nt)
               call write_log(message)
            endif
         endif
       enddo
      endif                     ! present(mtsum_init)

      end subroutine global_conservation

!----------------------------------------------------------------------

   subroutine glissade_add_smb(nx,         ny,         &
                               nlyr,       ntracer,    &
                               dt,                     &
                               thck_layer, tracer,     &
                               acab,       tsfc,       &
                               bmlt,       tbed)
      ! Input/output arguments

      integer, intent(in) ::   &
         nx, ny,               &! horizontal array size
         nlyr,                 &! number of vertical layers
         ntracer

      real(dp), intent(in) ::   &
         dt                     ! time step (s)

      real(dp), dimension (nx,ny,nlyr), intent(inout) ::     &
         thck_layer             ! ice layer thickness

      real(dp), dimension (nx,ny,ntracer,nlyr), intent(inout) ::     &
         tracer                 ! tracer values

      real(dp), intent(in), dimension(nx,ny) :: &
         acab                   ! surface mass balance (m/s)

      real(dp), intent(in), dimension(nx,ny), optional :: &
         tsfc                   ! surface temperature (deg C)

      real(dp), intent(in), dimension(nx,ny) :: &
         bmlt                   ! basal melt rate (m/s)
                                ! > 0 for melting, < 0 for freeze-on

      real(dp), intent(in), dimension(nx,ny), optional :: &
         tbed                   ! basal temperature (deg C)

      ! Local variables

      real(dp), dimension(nx,ny,ntracer,nlyr) ::  &
         thck_tracer       ! thck_layer * tracer

      real(dp) :: sfc_accum, sfc_ablat  ! surface accumulation/ablation, from acab
      real(dp) :: bed_accum, bed_ablat  ! bed accumulation/ablation, from bmlt
      real(dp) :: new_bmlt              ! bmlt adjusted for residual surface ablation

      integer :: i, j, k

      character(len=100) :: message

!WHL - debug
      real(dp) :: init_thck, new_thck

      do j = 1, ny
         do i = 1, nx

            ! initialize accumulation/ablation terms
            sfc_accum = 0.d0
            sfc_ablat = 0.d0
            bed_accum = 0.d0
            bed_ablat = 0.d0
            
!WHL - debug
!            init_thck = sum(thck_layer(i,j,:))
!            if (i==idiag .and. j==jdiag) then
!               print*, 'starting thck 1 =', thck_layer(i,j,1)
!               print*, 'starting thck =', sum(thck_layer(i,j,:))
!               print*, 'acab*dt =', acab(i,j)*dt
!               print*, 'bmlt*dt =', bmlt(i,j)*dt
!               print*, 'thck + acab*dt - bmlt*dt =', sum(thck_layer(i,j,:)) + acab(i,j)*dt - bmlt(i,j)*dt
!            endif

            ! Add surface accumulation/ablation to ice thickness
            ! Also modify tracers conservatively.
            ! Assume tracer(:,:,1,:) is temperature and is always present
            ! Asumme tracer(:,:,2,:) is ice age and is optionally present

            if (acab(i,j) > 0.d0) then       ! accumulation, added to layer 1

               sfc_accum = acab(i,j)*dt

               ! temperature of top layer
               if (ntracer >= 1) then
                  thck_tracer(i,j,1,1) = thck_layer(i,j,1) * tracer(i,j,1,1)  &
                                       + sfc_accum * tsfc(i,j)
               endif

               ! ice age (= 0 for new accumulation)
               if (ntracer >= 2) then  
                  thck_tracer(i,j,2,1) = thck_layer(i,j,1) * tracer(i,j,2,1) 
                                     ! + sfc_accum * 0.d0

                  !TODO - Add other tracers here, as needed

               endif

               ! new top layer thickess
               thck_layer(i,j,1) = thck_layer(i,j,1) + sfc_accum

               ! new tracer values in top layer
               tracer(i,j,:,1) = thck_tracer(i,j,:,1) / thck_layer(i,j,1)

            elseif (acab(i,j) < 0.d0) then   ! ablation             
                                             ! allow for ablation of multiple layers

               !TODO - Consider a 2nd-order accurate scheme allowing for linear
               !       tracer gradients within layers.

               ! reduce ice thickness (tracer values will not change)

               sfc_ablat = -acab(i,j)*dt   ! positive by definition

               !WHL - debug
!               if (init_thck > 0.d0 .and. init_thck < sfc_ablat) then
!                  print*, 'Surface ablation of entire ice thickness for cell', i, j
!                  print*, 'thck, sfc_ablat:', init_thck, sfc_ablat
!               endif

               do k = 1, nlyr
                  if (sfc_ablat > thck_layer(i,j,k)) then
                     sfc_ablat = sfc_ablat - thck_layer(i,j,k)
                     thck_layer(i,j,k) = 0.d0
                     tracer(i,j,:,k) = 0.d0
                  else
                     thck_layer(i,j,k) = thck_layer(i,j,k) - sfc_ablat
                     sfc_ablat = 0.d0
                     exit
                  endif
               enddo
  
            endif  ! acab > 0

            !TODO - Figure out how to handle energy conservation if we are left with sfc_ablat > 0.
            !       Include in the heat flux passed back to CLM?
            !       (This doesn't matter for EISMINT-type runs where acab < 0 in many ice-free cells.
            !        But in coupled runs, acab < 0 implies that we should melt a certain mass of ice.)

            !WHL - debug
!            if (init_thck > 0.d0 .and. sfc_ablat > 0.d0) then
!               print*, 'Residual energy for sfc ablation: i, j, sfc_ablat =', i, j, sfc_ablat
!            endif
!            new_thck = sum(thck_layer(i,j,:))

            ! Note: It is theoretically possible that we could have residual energy remaining for surface
            ! ablation while ice is freezing on at the bed, in which case the surface ablation should
            ! be subtracted from the bed accumulation.  But I'm going to ignore this possibility.

            if (bmlt(i,j) < 0.d0) then       ! freeze-on, added to lowest layer

               bed_accum = -bmlt(i,j)*dt

               ! temperature of bottom layer
               if (ntracer >= 1) then
                  thck_tracer(i,j,1,nlyr) = thck_layer(i,j,nlyr) * tracer(i,j,1,nlyr)  &
                                          + bed_accum * tbed(i,j)
               endif

               ! ice age (= 0 for new accumulation)
               if (ntracer >= 2) then  
                  thck_tracer(i,j,2,nlyr) = thck_layer(i,j,nlyr) * tracer(i,j,2,nlyr) 
                                        ! + bed_accum * 0.d0

                  !TODO - Add other tracers here, as needed
               endif

               ! new bottom layer thickess
               thck_layer(i,j,nlyr) = thck_layer(i,j,nlyr) + bed_accum

               ! new tracer values in bottom layer
               tracer(i,j,:,nlyr) = thck_tracer(i,j,:,nlyr) / thck_layer(i,j,nlyr)

            elseif (bmlt(i,j) > 0.d0) then   ! basal melting             
                                             ! allow for melting of multiple layers

               !TODO - Consider a 2nd-order accurate scheme allowing for linear
               !       tracer gradients within layers.

               ! reduce ice thickness (tracer values will not change)

               bed_ablat = bmlt(i,j)*dt   ! positive by definition

               !WHL - debug
!               if (new_thck > 0.d0 .and. new_thck < bed_ablat) then
!                  print*, 'Basal ablation of entire ice thickness for cell', i, j
!                  print*, 'thck, bed_ablat =', new_thck, bed_ablat
!               endif

               do k = nlyr, 1, -1
                  if (bed_ablat > thck_layer(i,j,k)) then
                     bed_ablat = bed_ablat - thck_layer(i,j,k)
                     thck_layer(i,j,k) = 0.d0
                     tracer(i,j,:,k) = 0.d0
                  else
                     thck_layer(i,j,k) = thck_layer(i,j,k) - bed_ablat
                     bed_ablat = 0.d0
                     exit
                  endif
               enddo
  
               !TODO - Figure out how to handle energy conservation if we are left with bed_ablat > 0.
               !       Include in the heat flux passed back to CLM?

               !WHL - debug
!               if (new_thck > 0.d0 .and. bed_ablat > 0.d0) then
!                  print*, 'Residual basal melting potential: i, j, bed_ablat =', i, j, bed_ablat
!               endif

            endif  ! bmlt < 0

         enddo   ! nx
      enddo      ! ny

      end subroutine glissade_add_smb

!----------------------------------------------------------------------

    subroutine glissade_vertical_remap(nx,       ny,        &
                                       nlyr,     ntracer,   &
                                       sigma,    hlyr,      &
                                       trcr)
 
    ! Conservative remapping of tracer fields from one set of vertical 
    ! coordinates to another.  The remapping is first-order accurate.
    !
    ! TODO - Add a 2nd-order accurate vertical remapping scheme.
    !
    ! Author: William Lipscomb, LANL

    implicit none
 
    ! in-out arguments
 
    integer, intent(in) ::  &
         nx, ny,     &! number of cells in EW and NS directions
         nlyr,       &! number of vertical layers
         ntracer      ! number of tracer fields

    real(dp), dimension (nx, ny, nlyr), intent(inout) ::  &
         hlyr         ! layer thickness

    real(dp), dimension (nlyr+1), intent(in) ::  &
         sigma        ! sigma vertical coordinate (at layer interfaces)

    real(dp), dimension (nx, ny, ntracer, nlyr), intent(inout) ::   &
         trcr         ! tracer field to be remapped
                      ! tracer(k) = value at midpoint of layer k
 
    ! local variables
 
    integer :: i, j, k, k1, k2, nt
 
    real(dp), dimension (nx, ny, nlyr+1) ::  &
         zi1,        &! layer interfaces in old coordinate system
                      ! zi1(1) = 0. = value at top surface
                      ! zi1(k) = value at top of layer k
                      ! zi1(nlyr+1) = value at bottom surface (= 1 in sigma coordinates)
         zi2          ! layer interfaces in new coordinate system
                      ! Note: zi1(1) = zi2(1) = 0
                      !       zi1(nlyr+1) = zi2(nlyr+1) = 1

    real(dp), dimension(nx, ny) ::        &
         thck,       &! total thickness
         rthck        ! reciprocal of total thickness
 
    real(dp), dimension(nx,ny,ntracer,nlyr) ::       &
         htsum        ! sum of thickness*tracer in a layer         

    real(dp) :: zlo, zhi, hovlp

       !-----------------------------------------------------------------
       ! Compute total thickness and reciprocal thickness
       !-----------------------------------------------------------------

       thck(:,:) = 0.d0
       do k = 1, nlyr
          thck(:,:) = thck(:,:) + hlyr(:,:,k)
       enddo

       do j = 1, ny
          do i = 1, nx
             if (thck(i,j) > 0.d0) then
                rthck(i,j) = 1.d0/thck(i,j)
             else
                rthck(i,j) = 0.d0
             endif
          enddo
       enddo

       !-----------------------------------------------------------------
       ! Determine vertical coordinate zi1, given input layer thicknesses.
       ! These are the coordinates from which we start.
       !-----------------------------------------------------------------

       zi1(:,:,1) = 0.d0
       do k = 2, nlyr
          zi1(:,:,k) = zi1(:,:,k-1) +  hlyr(:,:,k-1)*rthck(:,:) 
       enddo
       zi1(:,:,nlyr+1) = 1.d0
                       
       !-----------------------------------------------------------------
       ! Set vertical coordinate zi2, given sigma.
       ! These are the coordinates to which we remap in the vertical.
       !-----------------------------------------------------------------

       zi2(:,:,1) = 0.d0
       do k = 2, nlyr
          zi2(:,:,k) = sigma(k)
       enddo
       zi2(:,:,nlyr+1) = 1.d0

       !-----------------------------------------------------------------
       ! Compute new layer thicknesses (zi2 coordinates)
       !-----------------------------------------------------------------

       do k = 1, nlyr
          hlyr(:,:,k) = (zi2(:,:,k+1) - zi2(:,:,k)) * thck(:,:)
       enddo

       !-----------------------------------------------------------------
       ! Compute sum of h*T for each new layer (k2) by integrating
       ! over the regions of overlap with old layers (k1).
       ! Note: It might be worth trying a more efficient
       !       search algorithm if the number of layers is large.
       !       This algorithm scales as nlyr^2.
       !       Also, may want to rearrange loop order if there are many tracers.
       !-----------------------------------------------------------------

       do k2 = 1, nlyr
          htsum(:,:,:,k2) = 0.d0 
          do k1 = 1, nlyr
             do nt = 1, ntracer
                do j = 1, ny
                   do i = 1, nx
                      zhi = min (zi1(i,j,k1+1), zi2(i,j,k2+1)) 
                      zlo = max (zi1(i,j,k1), zi2(i,j,k2))
                      hovlp = max (zhi-zlo, 0.d0) * thck(i,j)
                      htsum(i,j,nt,k2) = htsum(i,j,nt,k2)    &
                                       +  trcr(i,j,nt,k1) * hovlp
                   enddo   ! i
                enddo      ! j
             enddo         ! nt
          enddo            ! k1
       enddo               ! k2
 
       !-----------------------------------------------------------------
       ! Compute tracer values in new layers
       !-----------------------------------------------------------------
 
       do k = 1, nlyr
          do nt = 1, ntracer
             do j = 1, ny
                do i = 1, nx
                   if (hlyr(i,j,k) > 0.d0) then
                      trcr(i,j,nt,k) = htsum(i,j,nt,k) / hlyr(i,j,k)
                   else
                      trcr(i,j,nt,k) = 0.d0
                   endif
                enddo   ! i
             enddo      ! j
          enddo         ! nt
       enddo            ! k

    end subroutine glissade_vertical_remap

!=======================================================================
!BOP
!
! !IROUTINE: upwind_field - first-order upwind transport
!
! !INTERFACE:
!
      subroutine upwind_field (nx,       ny,         &
                               ilo, ihi, jlo, jhi,   &
                               dx,       dy,         &
                               dt,       phi,        &
                               uee,      vnn)
!
! !DESCRIPTION:
!
! upwind transport algorithm
!
! !REVISION HISTORY:
!
! Authors: Elizabeth Hunke and William Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer, intent (in) ::     &
         nx, ny             ,&! block dimensions
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real(dp), intent(in) ::         &
         dx, dy             ,&! x and y gridcell dimensions
         dt                   ! time step

      real(dp), dimension(nx,ny), &
         intent(inout) ::                       &
         phi                  ! scalar field

      real(dp), dimension(nx,ny),         &
         intent(in)::     &
         uee, vnn             ! cell edge velocities
!
!EOP
!
      integer ::     &
         i, j                   ! standard indices

      real(dp) ::        &
         upwind, y1, y2, a, h   ! function

      real(dp), dimension(nx,ny) ::  &
         worka, workb

    !-------------------------------------------------------------------
    ! Define upwind function
    !-------------------------------------------------------------------

      upwind(y1,y2,a,h) = 0.5d0*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

    !-------------------------------------------------------------------
    ! upwind transport
    !-------------------------------------------------------------------

      do j = jlo-1, jhi
      do i = ilo-1, ihi
         worka(i,j)= upwind(phi(i,j),phi(i+1,j),uee(i,j),dy)
         workb(i,j)= upwind(phi(i,j),phi(i,j+1),vnn(i,j),dx)
      enddo
      enddo

      do j = jlo, jhi
      do i = ilo, ihi
         phi(i,j) = phi(i,j) - ( worka(i,j)-worka(i-1,j)      &
                               + workb(i,j)-workb(i,j-1) )    &
                               / (dx*dy)
      enddo
      enddo

      end subroutine upwind_field

!=======================================================================

      end module glissade_transport

!=======================================================================
