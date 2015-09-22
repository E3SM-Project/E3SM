!=======================================================================
!BOP
!
! !MODULE: ice_dyn_evp - elastic-viscous-plastic sea ice dynamics model
!
! !DESCRIPTION:
!
! Elastic-viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. {\em J. Phys. Oceanogr.}, {\bf 27}, 1849--1867.
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. {\em Journal of Computational Physics}, {\bf 170},
! 18--38.
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere---Incorporation of Metric Terms. {\em Monthly Weather Review},
! {\bf 130}, 1848--1865.
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.
!
!
! !REVISION HISTORY:
!  SVN:$Id: ice_dyn_evp.F90 100 2008-01-29 00:25:32Z eclare $
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
! 2004: Block structure added by William Lipscomb
! 2005: Removed boundary calls for stress arrays (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)
! 
! !INTERFACE:
!
      module ice_dyn_evp
!
! !USES:
!
      use ice_kinds_mod
      use ice_fileunits
      use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
      use ice_domain_size
      use ice_constants
      use ice_exit, only: abort_ice
      use perf_mod,        only: t_startf, t_stopf, t_barrierf
!
!EOP
!
      implicit none
      save

      ! namelist parameters

      integer (kind=int_kind) :: &
         kdyn     , & ! type of dynamics ( 1 = evp )
         ndte         ! number of subcycles:  ndte=dt/dte
      logical (kind=log_kind) :: &
         maskhalo_dyn , &  ! turn on masked halo updates in subcycling
         maskhalo_stress , &  ! turn on masked halo updates in stress update for tripole
         splitcomm_dyn        ! turn on overlapping of halo update and computation in subcycling

      logical (kind=log_kind) :: &
         evp_damping  ! if true, use evp damping procedure

      ! other EVP parameters

      character (len=char_len) :: & 
         yield_curve  ! 'ellipse' ('teardrop' needs further testing)
                                                                      ! 
      real (kind=dbl_kind), parameter :: &
         dragw = dragio * rhow, &
                         ! drag coefficient for water on ice *rhow (kg/m^3)
         eyc = 0.36_dbl_kind, &
                         ! coefficient for calculating the parameter E
         cosw = c1   , & ! cos(ocean turning angle)  ! turning angle = 0
         sinw = c0   , & ! sin(ocean turning angle)  ! turning angle = 0
         a_min = p001, & ! minimum ice area
         m_min = p01     ! minimum ice mass (kg/m^2)

      real (kind=dbl_kind) :: &
         ecci     , & ! 1/e^2
         dtei     , & ! 1/dte, where dte is subcycling timestep (1/s)
         dte2T    , & ! dte/2T
         denom1   , & ! constants for stress equation
         denom2   , & !
         rcon         ! for damping criterion (kg/s)

      real (kind=dbl_kind), allocatable :: & 
         fcor_blk(:,:,:)   ! Coriolis parameter (1/s)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: evp - elastic-viscous-plastic dynamics driver
!
! !INTERFACE:
!
      subroutine evp (dt)
!
! !DESCRIPTION:
!
! Elastic-viscous-plastic dynamics driver
!
#ifdef CICE_IN_NEMO
! Wind stress is set during this routine from the values supplied
! via NEMO.  These values are supplied rotated on u grid and
! multiplied by aice.  strairxT = 0 in this case so operations in
! evp_prep1 are pointless but carried out to minimise code changes.
#endif
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_boundary
      use ice_blocks
      use ice_domain
      use ice_state
      use ice_flux
      use ice_grid
      use ice_timers
      use ice_mechred, only: ice_strength
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: & 
         ksub           , & ! subcycle step
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j

      integer (kind=int_kind), dimension(max_blocks) :: & 
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         tmass    , & ! total mass of ice and snow (kg/m^2)
         waterx   , & ! for ocean stress calculation, x (m/s)
         watery   , & ! for ocean stress calculation, y (m/s)
         forcex   , & ! work array: combined atm stress and ocn tilt, x
         forcey   , & ! work array: combined atm stress and ocn tilt, y
         aiu      , & ! ice fraction on u-grid
         umass    , & ! total mass of ice and snow (u grid)
         umassdtei    ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), allocatable :: fld2(:,:,:,:)

!      real (kind=dbl_kind), dimension(nx_block,ny_block,8,max_blocks):: &
!         str8         ! stress combinations for momentum equation

      real (kind=dbl_kind), allocatable ::  &
         str8(:,:,:,:)         ! stress combinations for momentum equation

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         icetmask   ! ice extent mask (T-cell)

      type (block) :: &
         this_block           ! block information for current block
      
      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         halomask     ! mask for masked halo creation
      type (ice_halo) :: &
         halo_info_mask          !  ghost cell update info

!      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

       ! This call is needed only if dt changes during runtime.
!echmod: automate this
!      call set_evp_parameters (dt)

      !-----------------------------------------------------------------
      ! boundary updates
      ! commented out because the ghost cells are freshly 
      ! updated after cleanup_itd
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (aice,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vice,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vsno,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_timer_stop(timer_bound)

!     call t_barrierf ('cice_dyn_evp_prep_BARRIER',MPI_COMM_ICE)

      call t_barrierf ('cice_evp_prep1_BARRIER',MPI_COMM_ICE)
      call t_startf   ('cice_evp_prep1')

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks

         do j = 1, ny_block 
         do i = 1, nx_block 
            rdg_conv (i,j,iblk) = c0 
            rdg_shear(i,j,iblk) = c0 
            divu (i,j,iblk) = c0 
            shear(i,j,iblk) = c0 
            prs_sig(i,j,iblk) = c0 
         enddo
         enddo

      !-----------------------------------------------------------------
      ! preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call evp_prep1 (nx_block,           ny_block,           & 
                         ilo, ihi,           jlo, jhi,           &
                         aice    (:,:,iblk), vice    (:,:,iblk), & 
                         vsno    (:,:,iblk), tmask   (:,:,iblk), & 
                         strairxT_accum(:,:,iblk), strairyT_accum(:,:,iblk), & 
                         strairx (:,:,iblk), strairy (:,:,iblk), & 
                         tmass   (:,:,iblk), icetmask(:,:,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call t_stopf   ('cice_evp_prep1')
      call t_barrierf('cice_evp_bound1_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_evp_bound1')

!      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (icetmask,          halo_info, &
                           field_loc_center,  field_type_scalar)
!      call ice_timer_stop(timer_bound)

      call t_stopf ('cice_evp_bound1')
      call t_barrierf ('cice_evp_convert_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_evp_convert')

      !-----------------------------------------------------------------
      ! convert fields from T to U grid
      !-----------------------------------------------------------------

      call to_ugrid(tmass,umass)
!EH and MM recommendation: 
      call to_ugrid(aice, aiu)        ! use this line when subcycling ice thermo
!      call to_ugrid(aice_init, aiu)  ! use this line without subcycling

#ifdef CICE_IN_NEMO
      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO
      ! This wind stress is rotated on u grid and multiplied by aice
      !----------------------------------------------------------------

      strairx(:,:,:) = strax(:,:,:)
      strairy(:,:,:) = stray(:,:,:)
#else
      call t2ugrid_vector(strairx)
      call t2ugrid_vector(strairy)
#endif

     call t_stopf ('cice_evp_convert')
     call t_barrierf ('cice_evp_prep2_BARRIER',MPI_COMM_ICE)
     call t_startf ('cice_evp_prep2')

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! more preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call evp_prep2 (nx_block,             ny_block,             & 
                         ilo, ihi,             jlo, jhi,             &
                         icellt(iblk),         icellu(iblk),         & 
                         indxti      (:,iblk), indxtj      (:,iblk), & 
                         indxui      (:,iblk), indxuj      (:,iblk), & 
                         aiu       (:,:,iblk), umass     (:,:,iblk), & 
                         umassdtei (:,:,iblk), fcor_blk  (:,:,iblk), & 
                         umask     (:,:,iblk),                       & 
                         uocn      (:,:,iblk), vocn      (:,:,iblk), & 
                         strairx   (:,:,iblk), strairy   (:,:,iblk), & 
                         ss_tltx   (:,:,iblk), ss_tlty   (:,:,iblk), &  
                         icetmask  (:,:,iblk), iceumask  (:,:,iblk), & 
                         fm        (:,:,iblk),                       & 
                         strtltx   (:,:,iblk), strtlty   (:,:,iblk), & 
                         strocnx   (:,:,iblk), strocny   (:,:,iblk), & 
                         strintx   (:,:,iblk), strinty   (:,:,iblk), & 
                         waterx    (:,:,iblk), watery    (:,:,iblk), & 
                         forcex    (:,:,iblk), forcey    (:,:,iblk), & 
                         stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                         stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                         stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                         uvel      (:,:,iblk), vvel      (:,:,iblk))

      !-----------------------------------------------------------------
      ! ice strength
      !-----------------------------------------------------------------

         call ice_strength (nx_block, ny_block,   & 
                            ilo, ihi, jlo, jhi,   &
                            icellt(iblk),         & 
                            indxti      (:,iblk), & 
                            indxtj      (:,iblk), & 
                            aice    (:,:,  iblk), & 
                            vice    (:,:,  iblk), & 
                            aice0   (:,:,  iblk), & 
                            aicen   (:,:,:,iblk), &  
                            vicen   (:,:,:,iblk), & 
                            strength(:,:,  iblk) )

      enddo  ! iblk
      !$OMP END PARALLEL DO

      call t_stopf ('cice_evp_prep2')
      call t_barrierf ('cice_evp_bound2_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_evp_bound2')

      allocate(fld2(nx_block,ny_block,2,max_blocks))

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1,nblocks
         fld2(1:nx_block,1:ny_block,1,iblk) = uvel(1:nx_block,1:ny_block,iblk)
         fld2(1:nx_block,1:ny_block,2,iblk) = vvel(1:nx_block,1:ny_block,iblk)
      enddo
      !$OMP END PARALLEL DO

      call ice_HaloUpdate (strength,           halo_info, &
                           field_loc_center,   field_type_scalar)
      ! velocities may have changed in evp_prep2
      call ice_HaloUpdate (fld2,               halo_info, &
                           field_loc_NEcorner, field_type_vector)

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1,nblocks
         uvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,1,iblk)
         vvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,2,iblk)
      enddo
      !$OMP END PARALLEL DO

     call t_stopf ('cice_evp_bound2')
     if (maskhalo_dyn) then
        call t_barrierf ('cice_evp_halocreate_BARRIER',MPI_COMM_ICE)
        call t_startf('cice_evp_halocreate')
        call ice_HaloMask(halo_info_mask, halo_info, icetmask)
        call t_stopf ('cice_evp_halocreate')
     endif

     allocate(str8(nx_block,ny_block,8,nblocks))

  if (splitcomm_dyn) then

      do ksub = 1,ndte        ! subcycling

         call t_barrierf ('cice_evp_subcycling_BARRIER',MPI_COMM_ICE)
         call t_startf ('cice_evp_subcycling')

         !-----------------------------------------------------------------
         ! send halo update, skip on first subcycle
         !-----------------------------------------------------------------

         if (ksub > 1) then
            call t_startf ('cice_evp_bound3i')
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1,nblocks
               fld2(1:nx_block,1:ny_block,1,iblk) = uvel(1:nx_block,1:ny_block,iblk)
               fld2(1:nx_block,1:ny_block,2,iblk) = vvel(1:nx_block,1:ny_block,iblk)
            enddo
            !$OMP END PARALLEL DO
            call t_stopf ('cice_evp_bound3i')
            call t_startf ('cice_evp_bound3s')
            if (maskhalo_dyn) then
               call ice_HaloUpdate (fld2,               halo_info_mask, &
                                    field_loc_NEcorner, field_type_vector, mode='send')
            else
               call ice_HaloUpdate (fld2,               halo_info, &
                                    field_loc_NEcorner, field_type_vector, mode='send')
            endif
            call t_stopf ('cice_evp_bound3s')
         endif

         !-----------------------------------------------------------------
         ! stress tensor equation, total surface stress, phase 1 = non edge pts
         !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            call t_startf ('cice_evp_stress1')
            call stress (1, ilo, ihi, jlo, jhi,                      &
                         nx_block,             ny_block,             & 
                         ksub,                 icellt(iblk),         & 
                         indxti      (:,iblk), indxtj      (:,iblk), & 
                         uvel      (:,:,iblk), vvel      (:,:,iblk), &     
                         dxt       (:,:,iblk), dyt       (:,:,iblk), & 
                         dxhy      (:,:,iblk), dyhx      (:,:,iblk), & 
                         cxp       (:,:,iblk), cyp       (:,:,iblk), & 
                         cxm       (:,:,iblk), cym       (:,:,iblk), & 
                         tarear    (:,:,iblk), tinyarea  (:,:,iblk), & 
                         strength  (:,:,iblk),                       & 
                         stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                         stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                         stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                         shear     (:,:,iblk), divu      (:,:,iblk), & 
                         prs_sig   (:,:,iblk),                       & 
                         rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), & 
                         str8    (:,:,:,iblk) )
            call t_stopf ('cice_evp_stress1')
         enddo
         !$OMP END PARALLEL DO

         !-----------------------------------------------------------------
         ! recv halo update, skip on first subcycle
         !-----------------------------------------------------------------

         if (ksub > 1) then
            call t_startf ('cice_evp_bound3r')
            if (maskhalo_dyn) then
               call ice_HaloUpdate (fld2,               halo_info_mask, &
                                    field_loc_NEcorner, field_type_vector, mode='recv')
            else
               call ice_HaloUpdate (fld2,               halo_info, &
                                    field_loc_NEcorner, field_type_vector, mode='recv')
            endif
            call t_stopf ('cice_evp_bound3r')

            call t_startf ('cice_evp_bound3c')
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1,nblocks
               uvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,1,iblk)
               vvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,2,iblk)
            enddo
            !$OMP END PARALLEL DO
            call t_stopf ('cice_evp_bound3c')
         endif

         !-----------------------------------------------------------------
         ! stress tensor equation, total surface stress, phase 2 = edge points
         !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            call t_startf ('cice_evp_stress2')
            call stress (2, ilo, ihi, jlo, jhi,                      &
                         nx_block,             ny_block,             & 
                         ksub,                 icellt(iblk),         & 
                         indxti      (:,iblk), indxtj      (:,iblk), & 
                         uvel      (:,:,iblk), vvel      (:,:,iblk), &     
                         dxt       (:,:,iblk), dyt       (:,:,iblk), & 
                         dxhy      (:,:,iblk), dyhx      (:,:,iblk), & 
                         cxp       (:,:,iblk), cyp       (:,:,iblk), & 
                         cxm       (:,:,iblk), cym       (:,:,iblk), & 
                         tarear    (:,:,iblk), tinyarea  (:,:,iblk), & 
                         strength  (:,:,iblk),                       & 
                         stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                         stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                         stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                         shear     (:,:,iblk), divu      (:,:,iblk), & 
                         prs_sig   (:,:,iblk),                       & 
                         rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), & 
                         str8    (:,:,:,iblk) )
            call t_stopf ('cice_evp_stress2')

         !-----------------------------------------------------------------
         ! momentum equation
         !-----------------------------------------------------------------

            call t_startf ('cice_evp_stepu')
            call stepu (nx_block,            ny_block,           & 
                        icellu       (iblk),                     & 
                        indxui     (:,iblk), indxuj    (:,iblk), & 
                        aiu      (:,:,iblk), str8  (:,:,:,iblk), & 
                        uocn     (:,:,iblk), vocn    (:,:,iblk), &     
                        waterx   (:,:,iblk), watery  (:,:,iblk), & 
                        forcex   (:,:,iblk), forcey  (:,:,iblk), & 
                        umassdtei(:,:,iblk), fm      (:,:,iblk), & 
                        uarear   (:,:,iblk),                     & 
                        strocnx  (:,:,iblk), strocny (:,:,iblk), & 
                        strintx  (:,:,iblk), strinty (:,:,iblk), & 
                        uvel     (:,:,iblk), vvel    (:,:,iblk))

            call t_stopf ('cice_evp_stepu')
         enddo
         !$OMP END PARALLEL DO

         call t_stopf ('cice_evp_subcycling')
        
      enddo                     ! subcycling

      !-----------------------------------------------------------------
      ! one final halo update on the velocities
      !-----------------------------------------------------------------

      call t_barrierf ('cice_evp_bound4_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_evp_bound4')
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1,nblocks
         fld2(1:nx_block,1:ny_block,1,iblk) = uvel(1:nx_block,1:ny_block,iblk)
         fld2(1:nx_block,1:ny_block,2,iblk) = vvel(1:nx_block,1:ny_block,iblk)
      enddo
      !$OMP END PARALLEL DO
      if (maskhalo_dyn) then
         call ice_HaloUpdate (fld2,               halo_info_mask, &
                              field_loc_NEcorner, field_type_vector)
       else
         call ice_HaloUpdate (fld2,               halo_info, &
                              field_loc_NEcorner, field_type_vector)
      endif
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1,nblocks
         uvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,1,iblk)
         vvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,2,iblk)
      enddo
      !$OMP END PARALLEL DO
      call t_stopf ('cice_evp_bound4')

  else

      do ksub = 1,ndte        ! subcycling

         call t_barrierf ('cice_evp_subcycling_BARRIER',MPI_COMM_ICE)
         call t_startf ('cice_evp_subcycling')

         !-----------------------------------------------------------------
         ! stress tensor equation, total surface stress, all gridcells
         !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            call t_startf ('cice_evp_stress0')
            call stress (0, ilo, ihi, jlo, jhi,                      &
                         nx_block,             ny_block,             & 
                         ksub,                 icellt(iblk),         & 
                         indxti      (:,iblk), indxtj      (:,iblk), & 
                         uvel      (:,:,iblk), vvel      (:,:,iblk), &     
                         dxt       (:,:,iblk), dyt       (:,:,iblk), & 
                         dxhy      (:,:,iblk), dyhx      (:,:,iblk), & 
                         cxp       (:,:,iblk), cyp       (:,:,iblk), & 
                         cxm       (:,:,iblk), cym       (:,:,iblk), & 
                         tarear    (:,:,iblk), tinyarea  (:,:,iblk), & 
                         strength  (:,:,iblk),                       & 
                         stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                         stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                         stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                         shear     (:,:,iblk), divu      (:,:,iblk), & 
                         prs_sig   (:,:,iblk),                       & 
                         rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), & 
                         str8    (:,:,:,iblk) )
            call t_stopf ('cice_evp_stress0')

         !-----------------------------------------------------------------
         ! momentum equation
         !-----------------------------------------------------------------

            call t_startf ('cice_evp_stepu')
            call stepu (nx_block,            ny_block,           & 
                        icellu       (iblk),                     & 
                        indxui     (:,iblk), indxuj    (:,iblk), & 
                        aiu      (:,:,iblk), str8  (:,:,:,iblk), & 
                        uocn     (:,:,iblk), vocn    (:,:,iblk), &     
                        waterx   (:,:,iblk), watery  (:,:,iblk), & 
                        forcex   (:,:,iblk), forcey  (:,:,iblk), & 
                        umassdtei(:,:,iblk), fm      (:,:,iblk), & 
                        uarear   (:,:,iblk),                     & 
                        strocnx  (:,:,iblk), strocny (:,:,iblk), & 
                        strintx  (:,:,iblk), strinty (:,:,iblk), & 
                        uvel     (:,:,iblk), vvel    (:,:,iblk))

            call t_stopf ('cice_evp_stepu')
         enddo
         !$OMP END PARALLEL DO

         !-----------------------------------------------------------------
         ! halo update
         !-----------------------------------------------------------------

            call t_startf ('cice_evp_bound3i')
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1,nblocks
               fld2(1:nx_block,1:ny_block,1,iblk) = uvel(1:nx_block,1:ny_block,iblk)
               fld2(1:nx_block,1:ny_block,2,iblk) = vvel(1:nx_block,1:ny_block,iblk)
            enddo
            !$OMP END PARALLEL DO
            call t_stopf ('cice_evp_bound3i')
            call t_startf ('cice_evp_bound3sr')
            if (maskhalo_dyn) then
               call ice_HaloUpdate (fld2,               halo_info_mask, &
                                    field_loc_NEcorner, field_type_vector)
            else
               call ice_HaloUpdate (fld2,               halo_info, &
                                    field_loc_NEcorner, field_type_vector)
            endif
            call t_stopf ('cice_evp_bound3sr')

            call t_startf ('cice_evp_bound3c')
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1,nblocks
               uvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,1,iblk)
               vvel(1:nx_block,1:ny_block,iblk) = fld2(1:nx_block,1:ny_block,2,iblk)
            enddo
            !$OMP END PARALLEL DO
            call t_stopf ('cice_evp_bound3c')

         call t_stopf ('cice_evp_subcycling')
        
      enddo                     ! subcycling

  endif   ! splitcomm subcycling

      deallocate(str8)
      deallocate(fld2)

      if (maskhalo_dyn) then
         call t_barrierf ('cice_evp_halodestroy_BARRIER',MPI_COMM_ICE)
         call t_startf ('cice_evp_halodestroy')
         call ice_HaloDestroy(halo_info_mask)
         call t_stopf ('cice_evp_halodestroy')
      endif

      call t_barrierf ('cice_evp_tpupd_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_evp_tpupd')

      ! Force symmetry across the tripole seam
      if (trim(grid_type) == 'tripole') then
      if (maskhalo_stress) then
         !-------------------------------------------------------
         ! set halomask to zero because ice_HaloMask always keeps
         ! local copies AND tripole zipper communication
         !-------------------------------------------------------
         call t_barrierf ('cice_evp_tpupdhc_BARRIER',MPI_COMM_ICE)
         call t_startf('cice_evp_tpupdhc')
         halomask = 0
         call ice_HaloMask(halo_info_mask, halo_info, halomask)
         call t_stopf ('cice_evp_tpupdhc')

         call t_barrierf ('cice_evp_tpupdhu_BARRIER',MPI_COMM_ICE)
         call t_startf('cice_evp_tpupdhu')
         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call t_stopf('cice_evp_tpupdhu')

         call t_barrierf ('cice_evp_tpupdhd_BARRIER',MPI_COMM_ICE)
         call t_startf ('cice_evp_tpupdhd')
         call ice_HaloDestroy(halo_info_mask)
         call t_stopf ('cice_evp_tpupdhd')
      else

         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info, &
                              field_loc_center,  field_type_scalar)

      endif
      endif

      call t_stopf ('cice_evp_tpupd')
      call t_barrierf ('cice_evp_finish_BARRIER',MPI_COMM_ICE)
      call t_startf ('cice_evp_finish')

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call evp_finish                               & 
              (nx_block,           ny_block,           & 
               icellu      (iblk),                     & 
               indxui    (:,iblk), indxuj    (:,iblk), & 
               uvel    (:,:,iblk), vvel    (:,:,iblk), & 
               uocn    (:,:,iblk), vocn    (:,:,iblk), & 
               aiu     (:,:,iblk),                     &
               strocnx (:,:,iblk), strocny (:,:,iblk), & 
               strocnxT(:,:,iblk), strocnyT(:,:,iblk))

      enddo
      !$OMP END PARALLEL DO

      call u2tgrid_vector(strocnxT)    ! shift
      call u2tgrid_vector(strocnyT)

!      call ice_timer_stop(timer_dynamics)    ! dynamics
      call t_stopf ('cice_evp_finish')

      end subroutine evp

!=======================================================================
!BOP
!
! !IROUTINE: init_evp - initialize parameters needed for evp dynamics
!
! !INTERFACE:
!
      subroutine init_evp (dt)
!
! !DESCRIPTION:
!
! Initialize parameters and variables needed for the evp dynamics
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_boundary
      use ice_blocks
      use ice_domain
      use ice_state
      use ice_flux
      use ice_grid
      use ice_fileunits
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, k, &
         iblk            ! block index

      real (kind=dbl_kind) :: &
         dte         , & ! subcycling timestep for EVP dynamics, s
         ecc         , & ! (ratio of major to minor ellipse axes)^2
         tdamp2          ! 2(wave damping time scale T)

      call set_evp_parameters (dt)

      if (my_task == master_task) then
         write(nu_diag,*) 'dt_dyn  = ',dt
         write(nu_diag,*) 'dte     = ',dt/real(ndte,kind=dbl_kind)
         write(nu_diag,*) 'tdamp   =', eyc*dt
      endif

      allocate(fcor_blk(nx_block,ny_block,max_blocks))

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block

         ! velocity
         uvel(i,j,iblk) = c0    ! m/s
         vvel(i,j,iblk) = c0    ! m/s

         ! strain rates
         divu (i,j,iblk) = c0
         shear(i,j,iblk) = c0
         rdg_conv (i,j,iblk) = c0
         rdg_shear(i,j,iblk) = c0

         ! Coriolis parameter
!!         fcor_blk(i,j,iblk) = 1.46e-4_dbl_kind ! Hibler 1979, N. Hem; 1/s
         fcor_blk(i,j,iblk) = c2*omega*sin(ULAT(i,j,iblk)) ! 1/s

         ! stress tensor,  kg/s^2
         stressp_1 (i,j,iblk) = c0
         stressp_2 (i,j,iblk) = c0
         stressp_3 (i,j,iblk) = c0
         stressp_4 (i,j,iblk) = c0
         stressm_1 (i,j,iblk) = c0
         stressm_2 (i,j,iblk) = c0
         stressm_3 (i,j,iblk) = c0
         stressm_4 (i,j,iblk) = c0
         stress12_1(i,j,iblk) = c0
         stress12_2(i,j,iblk) = c0
         stress12_3(i,j,iblk) = c0
         stress12_4(i,j,iblk) = c0

         ! ice extent mask on velocity points
         iceumask(i,j,iblk) = .false.

      enddo                     ! i
      enddo                     ! j
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      end subroutine init_evp

!=======================================================================
!BOP
!
! !IROUTINE: set_evp_parameters - set parameters for evp dynamics
!
! !INTERFACE:
!
      subroutine set_evp_parameters (dt)
!
! !DESCRIPTION:
!
! Set parameters needed for the evp dynamics.
! Note: This subroutine is currently called only during initialization.
!       If the dynamics time step can vary during runtime, it should
!        be called whenever the time step changes.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      real (kind=dbl_kind) :: &
         dte         , & ! subcycling timestep for EVP dynamics, s
         ecc         , & ! (ratio of major to minor ellipse axes)^2
         tdamp2          ! 2*(wave damping time scale T)

      ! elastic time step
      dte = dt/real(ndte,kind=dbl_kind)        ! s
      dtei = c1/dte              ! 1/s

      ! major/minor axis length ratio, squared
      ecc  = c4
      ecci = p25                  ! 1/ecc

      ! constants for stress equation
      tdamp2 = c2*eyc*dt                    ! s
      dte2T = dte/tdamp2                    ! ellipse (unitless)
      denom1 = c1/(c1+dte2T)
      denom2 = c1/(c1+dte2T*ecc)
      rcon = 1230._dbl_kind*eyc*dt*dtei**2  ! kg/s

      end subroutine set_evp_parameters

!=======================================================================
!BOP
!
! !IROUTINE: evp_prep1 - compute quantities needed for stress tensor and mom eqns
!
! !INTERFACE:
!
      subroutine evp_prep1 (nx_block,  ny_block, & 
                            ilo, ihi,  jlo, jhi, &
                            aice,      vice,     & 
                            vsno,      tmask,    & 
                            strairxT,  strairyT, & 
                            strairx,   strairy,  & 
                            tmass,     icetmask)

! !DESCRIPTION:
!
! Computes quantities needed in the stress tensor (sigma)
! and momentum (u) equations, but which do not change during
! the thermodynamics/transport time step:
! ice mass and ice extent masks
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(in) :: &
         aice    , & ! concentration of ice
         vice    , & ! volume per unit area of ice          (m)
         vsno    , & ! volume per unit area of snow         (m)
         strairxT, & ! stress on ice by air, x-direction
         strairyT    ! stress on ice by air, y-direction

      logical (kind=log_kind), dimension (nx_block,ny_block), & 
         intent(in) :: &
         tmask       ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(out) :: &
         strairx , & ! stress on ice by air, x-direction
         strairy , & ! stress on ice by air, y-direction
         tmass       ! total mass of ice and snow (kg/m^2)

      integer (kind=int_kind), dimension (nx_block,ny_block), & 
         intent(out) :: &
         icetmask    ! ice extent mask (T-cell)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j

      logical (kind=log_kind), dimension(nx_block,ny_block) :: &
         tmphm               ! temporary mask

      do j = 1, ny_block
      do i = 1, nx_block

      !-----------------------------------------------------------------
      ! total mass of ice and snow, centered in T-cell
      ! NOTE: vice and vsno must be up to date in all grid cells,
      !       including ghost cells
      !-----------------------------------------------------------------
         if (tmask(i,j)) then
            tmass(i,j) = (rhoi*vice(i,j) + rhos*vsno(i,j)) ! kg/m^2
         else
            tmass(i,j) = c0
         endif

      !-----------------------------------------------------------------
      ! ice extent mask (T-cells)
      !-----------------------------------------------------------------
         tmphm(i,j) = tmask(i,j) .and. (aice (i,j) > a_min) & 
                                 .and. (tmass(i,j) > m_min)

      !-----------------------------------------------------------------
      ! prep to convert to U grid
      !-----------------------------------------------------------------
         ! these quantities include the factor of aice needed for
         ! correct treatment of free drift
         strairx(i,j) = strairxT(i,j)
         strairy(i,j) = strairyT(i,j)

      !-----------------------------------------------------------------
      ! augmented mask (land + open ocean)
      !-----------------------------------------------------------------
         icetmask (i,j) = 0

      enddo
      enddo

      do j = jlo, jhi
      do i = ilo, ihi

         ! extend ice extent mask (T-cells) to points around pack
         if (tmphm(i-1,j+1) .or. tmphm(i,j+1) .or. tmphm(i+1,j+1) .or. & 
             tmphm(i-1,j)   .or. tmphm(i,j)   .or. tmphm(i+1,j)   .or. & 
             tmphm(i-1,j-1) .or. tmphm(i,j-1) .or. tmphm(i+1,j-1) ) then
            icetmask(i,j) = 1
         endif

         if (.not.tmask(i,j)) icetmask(i,j) = 0

      enddo
      enddo

      end subroutine evp_prep1

!=======================================================================
!BOP
!
! !IROUTINE: evp_prep2 - compute quantities needed for stress tensor and mom eqns
!
! !INTERFACE:
!
      subroutine evp_prep2 (nx_block,   ny_block,   & 
                            ilo, ihi,  jlo, jhi,    &
                            icellt,     icellu,     & 
                            indxti,     indxtj,     & 
                            indxui,     indxuj,     & 
                            aiu,        umass,      & 
                            umassdtei,  fcor,       & 
                            umask,                  & 
                            uocn,       vocn,       & 
                            strairx,    strairy,    & 
                            ss_tltx,    ss_tlty,    &  
                            icetmask,   iceumask,   & 
                            fm,                     & 
                            strtltx,    strtlty,    & 
                            strocnx,    strocny,    &
                            strintx,    strinty,    &
                            waterx,     watery,     & 
                            forcex,     forcey,     &     
                            stressp_1,  stressp_2,  &   
                            stressp_3,  stressp_4,  & 
                            stressm_1,  stressm_2,  & 
                            stressm_3,  stressm_4,  & 
                            stress12_1, stress12_2, & 
                            stress12_3, stress12_4, & 
                            uvel,       vvel)

! !DESCRIPTION:
!
! Computes quantities needed in the stress tensor (sigma)
! and momentum (u) equations, but which do not change during
! the thermodynamics/transport time step:
! --wind stress shift to U grid,
! --ice mass and ice extent masks,
! initializes ice velocity for new points to ocean sfc current
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      integer (kind=int_kind), intent(out) :: &
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(out) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      logical (kind=log_kind), dimension (nx_block,ny_block), & 
         intent(in) :: &
         umask       ! land/boundary mask, thickness (U-cell)

      integer (kind=int_kind), dimension (nx_block,ny_block), & 
         intent(in) :: &
         icetmask    ! ice extent mask (T-cell)

      logical (kind=log_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         iceumask    ! ice extent mask (U-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aiu     , & ! ice fraction on u-grid
         umass   , & ! total mass of ice and snow (u grid)
         fcor    , & ! Coriolis parameter (1/s)
         strairx , & ! stress on ice by air, x-direction
         strairy , & ! stress on ice by air, y-direction
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         ss_tltx , & ! sea surface slope, x-direction (m/m)
         ss_tlty     ! sea surface slope, y-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(out) :: &
         umassdtei,& ! mass of U-cell/dte (kg/m^2 s)
         waterx  , & ! for ocean stress calculation, x (m/s)
         watery  , & ! for ocean stress calculation, y (m/s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         forcey      ! work array: combined atm stress and ocn tilt, y

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         fm      , & ! Coriolis param. * mass in U-cell (kg/s)
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4, & ! sigma12
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         strtltx , & ! stress due to sea surface slope, x-direction
         strtlty , & ! stress due to sea surface slope, y-direction
         strocnx , & ! ice-ocean stress, x-direction
         strocny , & ! ice-ocean stress, y-direction
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         strinty     ! divergence of internal ice stress, y (N/m^2)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij

      logical (kind=log_kind), dimension(nx_block,ny_block) :: &
         iceumask_old      ! old-time iceumask

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         waterx   (i,j) = c0
         watery   (i,j) = c0
         forcex   (i,j) = c0
         forcey   (i,j) = c0
         umassdtei(i,j) = c0

         if (icetmask(i,j)==0) then
            stressp_1 (i,j) = c0
            stressp_2 (i,j) = c0
            stressp_3 (i,j) = c0
            stressp_4 (i,j) = c0
            stressm_1 (i,j) = c0
            stressm_2 (i,j) = c0
            stressm_3 (i,j) = c0
            stressm_4 (i,j) = c0
            stress12_1(i,j) = c0
            stress12_2(i,j) = c0
            stress12_3(i,j) = c0
            stress12_4(i,j) = c0
         endif                  ! icetmask
      enddo                     ! i
      enddo                     ! j

      !-----------------------------------------------------------------
      ! Identify cells where icetmask = 1
      ! Note: The icellt mask includes north and east ghost cells
      !       where stresses are needed.
      !-----------------------------------------------------------------

      icellt = 0
      do j = jlo, jhi+1
      do i = ilo, ihi+1
         if (icetmask(i,j) == 1) then
            icellt = icellt + 1
            indxti(icellt) = i
            indxtj(icellt) = j
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Define iceumask
      ! Identify cells where iceumask is true
      ! Initialize velocity where needed
      !-----------------------------------------------------------------

      icellu = 0
      do j = jlo, jhi
      do i = ilo, ihi

         ! ice extent mask (U-cells)
         iceumask_old(i,j) = iceumask(i,j) ! save
         iceumask(i,j) = (umask(i,j)) .and. (aiu  (i,j) > a_min) & 
                                      .and. (umass(i,j) > m_min)

         if (iceumask(i,j)) then
            icellu = icellu + 1
            indxui(icellu) = i
            indxuj(icellu) = j

            ! initialize velocity for new ice points to ocean sfc current
            if (.not. iceumask_old(i,j)) then
               uvel(i,j) = uocn(i,j)
               vvel(i,j) = vocn(i,j)
            endif
         else
            ! set velocity and stresses to zero for masked-out points
            uvel(i,j)    = c0
            vvel(i,j)    = c0
            strintx(i,j) = c0
            strinty(i,j) = c0
            strocnx(i,j) = c0
            strocny(i,j) = c0
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Define variables for momentum equation
      !-----------------------------------------------------------------

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         umassdtei(i,j) = umass(i,j)*dtei ! m/dte, kg/m^2 s
         fm(i,j) = fcor(i,j)*umass(i,j)   ! Coriolis * mass

         ! for ocean stress
         waterx(i,j) = uocn(i,j)*cosw - vocn(i,j)*sinw
         watery(i,j) = vocn(i,j)*cosw + uocn(i,j)*sinw

         ! combine tilt with wind stress
#ifndef coupled
         ! calculate tilt from geostrophic currents if needed
         strtltx(i,j) = -fm(i,j)*vocn(i,j)
         strtlty(i,j) =  fm(i,j)*uocn(i,j)
#else
         strtltx(i,j) = -gravit*umass(i,j)*ss_tltx(i,j)
         strtlty(i,j) = -gravit*umass(i,j)*ss_tlty(i,j)
#endif
         forcex(i,j) = strairx(i,j) + strtltx(i,j)
         forcey(i,j) = strairy(i,j) + strtlty(i,j)
      enddo

      end subroutine evp_prep2

!=======================================================================
!BOP
!
! !IROUTINE: stress - computes strain rates and internal stress components
!
! !INTERFACE:
!
      subroutine stress (phase, ilo, ihi, jlo, jhi, &
                         nx_block,   ny_block,   & 
                         ksub,       icellt,     & 
                         indxti,     indxtj,     & 
                         uvel,       vvel,       & 
                         dxt,        dyt,        & 
                         dxhy,       dyhx,       & 
                         cxp,        cyp,        & 
                         cxm,        cym,        & 
                         tarear,     tinyarea,   & 
                         strength,               & 
                         stressp_1,  stressp_2,  & 
                         stressp_3,  stressp_4,  & 
                         stressm_1,  stressm_2,  & 
                         stressm_3,  stressm_4,  & 
                         stress12_1, stress12_2, & 
                         stress12_3, stress12_4, & 
                         shear,      divu,       & 
                         prs_sig,                & 
                         rdg_conv,   rdg_shear,  & 
                         str8 )
!
! !DESCRIPTION:
!
! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: & 
         phase,              & ! phase
         ilo, ihi, jlo, jhi, & ! block dimensions
         nx_block, ny_block, & ! block dimensions
         ksub              , & ! subcycling step
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         strength , & ! ice strength (N/m)
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm      , & ! 0.5*HTN - 1.5*HTN
         tarear   , & ! 1/tarea
         tinyarea     ! puny*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         prs_sig  , & ! replacement pressure, for stress calc
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
         intent(out) :: &
         str8         ! stress combinations
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij, chunk

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
        c0ne, c0nw, c0se, c0sw                    , & ! useful combinations
        c1ne, c1nw, c1se, c1sw                    , &
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp, tmp

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      if (phase == 0 .or. phase == 1 .or. phase == 2) then
         ! valid, 0=all cells, 1=interior cells, 2=edge cells
      else
         call abort_ice ('ice_dyn stress: illegal phase')
      endif

      if (phase == 1) then
         str8(:,:,:) = c0
      endif

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

!     if ((phase == 1 .and. (i >  ilo .and. i <  ihi) .and. (j >  jlo .and. j <  jhi)) .or. &
!         (phase == 2 .and. (i <= ilo .or.  i >= ihi) .and. (j <= jlo .or.  j >= jhi))) then

     if (i >  ilo .and. i <  ihi .and. j > jlo .and. j < jhi) then
        chunk = 1
     else
        chunk = 2
     endif

     if (phase == 0 .or. phase == chunk) then
      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         ! divergence  =  e_11 + e_22
         divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                   + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
         divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                   + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
         divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                   + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
         divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                   + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )

         ! tension strain rate  =  e_11 - e_22
         tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                   +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
         tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                   +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
         tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                   +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
         tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                   +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )

         ! shearing strain rate  =  e_12
         shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  ) &
                 -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
         shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  ) &
                 -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
         shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1) &
                 -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
         shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1) &
                 -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )

         ! Delta (in the denominator of zeta, eta)
         Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
         Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
         Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
         Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))

      !-----------------------------------------------------------------
      ! on last subcycle, save quantities for mechanical redistribution
      !-----------------------------------------------------------------
         if (ksub == ndte) then
            divu(i,j) = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
            tmp = p25*(Deltane + Deltanw + Deltase + Deltasw)   * tarear(i,j)
            rdg_conv(i,j)  = -min(divu(i,j),c0)
            rdg_shear(i,j) = p5*(tmp-abs(divu(i,j))) 

            ! diagnostic only
            ! shear = sqrt(tension**2 + shearing**2)
            shear(i,j) = p25*tarear(i,j)*sqrt( &
                 (tensionne + tensionnw + tensionse + tensionsw)**2 &
                +  (shearne +   shearnw +   shearse +   shearsw)**2)

         endif

      !-----------------------------------------------------------------
      ! replacement pressure/Delta                   ! kg/s
      ! save replacement pressure for principal stress calculation
      !-----------------------------------------------------------------
         if (evp_damping) then
            ! enforce damping criterion
            c0ne = min(strength(i,j)/max(Deltane,c4*tinyarea(i,j)),rcon)
            c0nw = min(strength(i,j)/max(Deltanw,c4*tinyarea(i,j)),rcon)
            c0sw = min(strength(i,j)/max(Deltasw,c4*tinyarea(i,j)),rcon)
            c0se = min(strength(i,j)/max(Deltase,c4*tinyarea(i,j)),rcon)
            prs_sig(i,j) = strength(i,j)* &
                           Deltane/max(Deltane,c4*tinyarea(i,j)) ! ne
         else
            ! original version
            c0ne = strength(i,j)/max(Deltane,tinyarea(i,j))
            c0nw = strength(i,j)/max(Deltanw,tinyarea(i,j))
            c0sw = strength(i,j)/max(Deltasw,tinyarea(i,j))
            c0se = strength(i,j)/max(Deltase,tinyarea(i,j))
            prs_sig(i,j) = c0ne*Deltane ! northeast
         endif

         c1ne = c0ne*dte2T
         c1nw = c0nw*dte2T
         c1sw = c0sw*dte2T
         c1se = c0se*dte2T

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

         stressp_1(i,j) = (stressp_1(i,j) + c1ne*(divune - Deltane)) &
                          * denom1
         stressp_2(i,j) = (stressp_2(i,j) + c1nw*(divunw - Deltanw)) &
                          * denom1
         stressp_3(i,j) = (stressp_3(i,j) + c1sw*(divusw - Deltasw)) &
                          * denom1
         stressp_4(i,j) = (stressp_4(i,j) + c1se*(divuse - Deltase)) &
                          * denom1

         stressm_1(i,j) = (stressm_1(i,j) + c1ne*tensionne) * denom2
         stressm_2(i,j) = (stressm_2(i,j) + c1nw*tensionnw) * denom2
         stressm_3(i,j) = (stressm_3(i,j) + c1sw*tensionsw) * denom2
         stressm_4(i,j) = (stressm_4(i,j) + c1se*tensionse) * denom2

         stress12_1(i,j) = (stress12_1(i,j) + c1ne*shearne*p5) * denom2
         stress12_2(i,j) = (stress12_2(i,j) + c1nw*shearnw*p5) * denom2
         stress12_3(i,j) = (stress12_3(i,j) + c1sw*shearsw*p5) * denom2
         stress12_4(i,j) = (stress12_4(i,j) + c1se*shearse*p5) * denom2

      !-----------------------------------------------------------------
      ! Eliminate underflows.
      ! The following code is commented out because it is relatively 
      ! expensive and most compilers include a flag that accomplishes
      ! the same thing more efficiently.  This code is cheaper than
      ! handling underflows if the compiler lacks a flag; uncomment
      ! it in that case.  The compiler flag is often described with the 
      ! phrase "flush to zero".
      !-----------------------------------------------------------------

!      stressp_1(i,j) = sign(max(abs(stressp_1(i,j)),puny),stressp_1(i,j))
!      stressp_2(i,j) = sign(max(abs(stressp_2(i,j)),puny),stressp_2(i,j))
!      stressp_3(i,j) = sign(max(abs(stressp_3(i,j)),puny),stressp_3(i,j))
!      stressp_4(i,j) = sign(max(abs(stressp_4(i,j)),puny),stressp_4(i,j))

!      stressm_1(i,j) = sign(max(abs(stressm_1(i,j)),puny),stressm_1(i,j))
!      stressm_2(i,j) = sign(max(abs(stressm_2(i,j)),puny),stressm_2(i,j))
!      stressm_3(i,j) = sign(max(abs(stressm_3(i,j)),puny),stressm_3(i,j))
!      stressm_4(i,j) = sign(max(abs(stressm_4(i,j)),puny),stressm_4(i,j))

!      stress12_1(i,j) = sign(max(abs(stress12_1(i,j)),puny),stress12_1(i,j))
!      stress12_2(i,j) = sign(max(abs(stress12_2(i,j)),puny),stress12_2(i,j))
!      stress12_3(i,j) = sign(max(abs(stress12_3(i,j)),puny),stress12_3(i,j))
!      stress12_4(i,j) = sign(max(abs(stress12_4(i,j)),puny),stress12_4(i,j))

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1(i,j) + stressp_2(i,j)
         ssigps  = stressp_3(i,j) + stressp_4(i,j)
         ssigpe  = stressp_1(i,j) + stressp_4(i,j)
         ssigpw  = stressp_2(i,j) + stressp_3(i,j)
         ssigp1  =(stressp_1(i,j) + stressp_3(i,j))*p055
         ssigp2  =(stressp_2(i,j) + stressp_4(i,j))*p055

         ssigmn  = stressm_1(i,j) + stressm_2(i,j)
         ssigms  = stressm_3(i,j) + stressm_4(i,j)
         ssigme  = stressm_1(i,j) + stressm_4(i,j)
         ssigmw  = stressm_2(i,j) + stressm_3(i,j)
         ssigm1  =(stressm_1(i,j) + stressm_3(i,j))*p055
         ssigm2  =(stressm_2(i,j) + stressm_4(i,j))*p055

         ssig12n = stress12_1(i,j) + stress12_2(i,j)
         ssig12s = stress12_3(i,j) + stress12_4(i,j)
         ssig12e = stress12_1(i,j) + stress12_4(i,j)
         ssig12w = stress12_2(i,j) + stress12_3(i,j)
         ssig121 =(stress12_1(i,j) + stress12_3(i,j))*p111
         ssig122 =(stress12_2(i,j) + stress12_4(i,j))*p111

         csigpne = p111*stressp_1(i,j) + ssigp2 + p027*stressp_3(i,j)
         csigpnw = p111*stressp_2(i,j) + ssigp1 + p027*stressp_4(i,j)
         csigpsw = p111*stressp_3(i,j) + ssigp2 + p027*stressp_1(i,j)
         csigpse = p111*stressp_4(i,j) + ssigp1 + p027*stressp_2(i,j)
         
         csigmne = p111*stressm_1(i,j) + ssigm2 + p027*stressm_3(i,j)
         csigmnw = p111*stressm_2(i,j) + ssigm1 + p027*stressm_4(i,j)
         csigmsw = p111*stressm_3(i,j) + ssigm2 + p027*stressm_1(i,j)
         csigmse = p111*stressm_4(i,j) + ssigm1 + p027*stressm_2(i,j)
         
         csig12ne = p222*stress12_1(i,j) + ssig122 &
                  + p055*stress12_3(i,j)
         csig12nw = p222*stress12_2(i,j) + ssig121 &
                  + p055*stress12_4(i,j)
         csig12sw = p222*stress12_3(i,j) + ssig122 &
                  + p055*stress12_1(i,j)
         csig12se = p222*stress12_4(i,j) + ssig121 &
                  + p055*stress12_2(i,j)

         str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         str8(i,j,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         str8(i,j,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str8(i,j,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         str8(i,j,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str8(i,j,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         str8(i,j,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str8(i,j,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         str8(i,j,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

      endif                     ! phase
      enddo                     ! ij

      end subroutine stress

!=======================================================================
!BOP
!
! !IROUTINE: stepu - integrates mom eqn for u,v
!
! !INTERFACE:
!
      subroutine stepu (nx_block,   ny_block, &
                        icellu,               &
                        indxui,     indxuj,   &
                        aiu,        str8,      &
                        uocn,       vocn,     &
                        waterx,     watery,   &
                        forcex,     forcey,   &
                        umassdtei,  fm,       &
                        uarear,               &
                        strocnx,    strocny,  &
                        strintx,    strinty,  &
                        uvel,       vvel)
!
! !DESCRIPTION:
!
! Calculation of the surface stresses
! Integration of the momentum equation to find velocity (u,v)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aiu     , & ! ice fraction on u-grid
         waterx  , & ! for ocean stress calculation, x (m/s)
         watery  , & ! for ocean stress calculation, y (m/s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         forcey  , & ! work array: combined atm stress and ocn tilt, y
         umassdtei,& ! mass of U-cell/dte (kg/m^2 s)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         fm      , & ! Coriolis param. * mass in U-cell (kg/s)
         uarear      ! 1/uarea

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), &
         intent(in) :: &
         str8

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel        ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         strocnx , & ! ice-ocean stress, x-direction
         strocny , & ! ice-ocean stress, y-direction
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         strinty     ! divergence of internal ice stress, y (N/m^2)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         uold, vold        , & ! old-time uvel, vvel
         vrel              , & ! relative ice-ocean velocity
         cca,ccb,ab2,cc1,cc2,& ! intermediate variables
         taux, tauy            ! part of ocean stress term          

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         uold = uvel(i,j)
         vold = vvel(i,j)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiu(i,j)*dragw*sqrt((uocn(i,j) - uold)**2 + &
                                    (vocn(i,j) - vold)**2)  ! m/s
         ! ice/ocean stress
         taux = vrel*waterx(i,j) ! NOTE this is not the entire
         tauy = vrel*watery(i,j) ! ocn stress term

         ! alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
         cca = umassdtei(i,j) + vrel * cosw      ! alpha, kg/m^2 s
         ccb = fm(i,j)        + vrel * sinw      ! beta,  kg/m^2 s
         ab2 = cca**2 + ccb**2

         ! divergence of the internal stress tensor
         strintx(i,j) = uarear(i,j)* &
             (str8(i,j,1) + str8(i+1,j,2) + str8(i,j+1,3) + str8(i+1,j+1,4))
         strinty(i,j) = uarear(i,j)* &
             (str8(i,j,5) + str8(i,j+1,6) + str8(i+1,j,7) + str8(i+1,j+1,8))

         ! finally, the velocity components
         cc1 = strintx(i,j) + forcex(i,j) + taux &
             + umassdtei(i,j)*uold
         cc2 = strinty(i,j) + forcey(i,j) + tauy &
             + umassdtei(i,j)*vold

         uvel(i,j) = (cca*cc1 + ccb*cc2) / ab2 ! m/s
         vvel(i,j) = (cca*cc2 - ccb*cc1) / ab2

      !-----------------------------------------------------------------
      ! ocean-ice stress for coupling
      ! here, strocn includes the factor of aice
      !-----------------------------------------------------------------
         strocnx(i,j) = taux
         strocny(i,j) = tauy

      enddo                     ! ij

      end subroutine stepu

!=======================================================================
!BOP
!
! !IROUTINE: evp_finish - calculates ice-ocean stress
!
! !INTERFACE:
!
      subroutine evp_finish (nx_block, ny_block, &
                             icellu,             &
                             indxui,   indxuj,   &
                             uvel,     vvel,     &
                             uocn,     vocn,     &
                             aiu,                &
                             strocnx,  strocny,  &
                             strocnxT, strocnyT) 
!
! !DESCRIPTION:
!
! Calculation of the ice-ocean stress.
! ...the sign will be reversed later...
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , &  ! ocean current, y-direction (m/s)
         aiu         ! ice fraction on u-grid

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         strocnx , & ! ice-ocean stress, x-direction
         strocny , & ! ice-ocean stress, y-direction
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT    ! ice-ocean stress, y-direction
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: vrel

      do j = 1, ny_block
      do i = 1, nx_block
         strocnxT(i,j) = c0
         strocnyT(i,j) = c0
      enddo
      enddo

      ! ocean-ice stress for coupling
      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

!maltrud implemeent Adrians fix
!        vrel = dragw*sqrt((uocn(i,j) - uvel(i,j))**2 + &
!                          (vocn(i,j) - vvel(i,j))**2)  ! m/s
!        strocnx(i,j) = strocnx(i,j) &
!                     - vrel*(uvel(i,j)*cosw - vvel(i,j)*sinw) * aiu(i,j)
!        strocny(i,j) = strocny(i,j) &
!                     - vrel*(vvel(i,j)*cosw + uvel(i,j)*sinw) * aiu(i,j)
         vrel = dragw*sqrt((uocn(i,j) - uvel(i,j))**2 + &
                           (vocn(i,j) - vvel(i,j))**2)  ! m/s
         strocnx(i,j) = vrel*(uocn(i,j)*cosw - vocn(i,j)*sinw) * aiu(i,j) &
                      - vrel*(uvel(i,j)*cosw - vvel(i,j)*sinw) * aiu(i,j)
         strocny(i,j) = vrel*(vocn(i,j)*cosw + uocn(i,j)*sinw) * aiu(i,j) &
                      - vrel*(vvel(i,j)*cosw + uvel(i,j)*sinw) * aiu(i,j)

         ! Prepare to convert to T grid
         ! divide by aice for coupling
         strocnxT(i,j) = strocnx(i,j) / aiu(i,j)
         strocnyT(i,j) = strocny(i,j) / aiu(i,j)
      enddo

      end subroutine evp_finish

!=======================================================================
!BOP
!
! !IROUTINE: principal_stress - computes principal stress for yield curve
!
! !INTERFACE:
!
      subroutine principal_stress(nx_block,   ny_block,  &
                                  stressp_1,  stressm_1, &
                                  stress12_1, prs_sig,   &
                                  sig1,       sig2)
!
! !DESCRIPTION:
!
! Computes principal stresses for comparison with the theoretical
! yield curve; northeast values
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block  ! block dimensions

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         stressp_1 , & ! sigma11 + sigma22
         stressm_1 , & ! sigma11 - sigma22
         stress12_1, & ! sigma12
         prs_sig       ! replacement pressure, for stress calc

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         sig1    , & ! principal stress component
         sig2        ! principal stress component
!
!EOP
!
      integer (kind=int_kind) :: i, j

      do j = 1, ny_block
      do i = 1, nx_block
         if (prs_sig(i,j) > puny) then
            sig1(i,j) = (p5*(stressp_1(i,j) &
                      + sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2))) &
                      / prs_sig(i,j)
            sig2(i,j) = (p5*(stressp_1(i,j) &
                      - sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2))) &
                      / prs_sig(i,j)
         else
            sig1(i,j) = spval_dbl
            sig2(i,j) = spval_dbl
         endif
      enddo
      enddo

      end subroutine principal_stress

!=======================================================================

      end module ice_dyn_evp

!=======================================================================
