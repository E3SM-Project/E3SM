!=======================================================================
!BOP
!
! !MODULE: ice_grid - spatial grids, masks and boundary conditions
!
! !DESCRIPTION:
!
! Spatial grids, masks, and boundary conditions
!
! !REVISION HISTORY:
!  SVN:$Id: ice_grid.F90 152 2008-09-24 20:48:59Z eclare $
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!          Tony Craig, NCAR
!
! 2004: Block structure added by William Lipscomb
!       init_grid split into two parts as in POP 2.0
!       Boundary update routines replaced by POP versions
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2007: Option to read from netcdf files (A. Keen, Met Office)
!       Grid reading routines reworked by E. Hunke for boundary values
!
! !INTERFACE:
!
      module ice_grid
!
! !USES:
!
      use ice_kinds_mod
      use ice_boundary
      use ice_communicate, only: my_task, master_task
      use ice_constants
      use ice_blocks
      use ice_domain_size
      use ice_domain
      use ice_fileunits
      use ice_gather_scatter
      use ice_read_write
      use ice_timers
      use ice_probability
      use ice_exit
!
!EOP
!
      implicit none
!echmod      save

      character (len=char_len_long), save :: &
         grid_format  , & ! file format ('bin'=binary or 'nc'=netcdf)
         grid_file    , & !  input file for POP grid info
         kmt_file     , & !  input file for POP grid info
         gridcpl_file , & !  input file for POP coupling grid info
         grid_type        !  current options are rectangular (default),
                          !  displaced_pole, tripole, panarctic

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         dxt    , & ! width of T-cell through the middle (m)
         dyt    , & ! height of T-cell through the middle (m)
         dxu    , & ! width of U-cell through the middle (m)
         dyu    , & ! height of U-cell through the middle (m)
         HTE    , & ! length of eastern edge of T-cell (m)
         HTN    , & ! length of northern edge of T-cell (m)
         tarea  , & ! area of T-cell (m^2)
         uarea  , & ! area of U-cell (m^2)
         tarear , & ! 1/tarea
         uarear , & ! 1/uarea
         tinyarea,& ! puny*tarea
         tarean , & ! area of NH T-cells
         tareas , & ! area of SH T-cells
         ULON   , & ! longitude of velocity pts (radians)
         ULAT   , & ! latitude of velocity pts (radians)
         TLON   , & ! longitude of temp pts (radians)
         TLAT   , & ! latitude of temp pts (radians)
         ANGLE  , & ! for conversions between POP grid and lat/lon
         ANGLET     ! ANGLE converted to T-cells

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         cyp    , & ! 1.5*HTE - 0.5*HTE
         cxp    , & ! 1.5*HTN - 0.5*HTN
         cym    , & ! 0.5*HTE - 1.5*HTE
         cxm    , & ! 0.5*HTN - 1.5*HTN
         dxhy   , & ! 0.5*(HTE - HTE)
         dyhx       ! 0.5*(HTN - HTN)

      ! Corners of grid boxes for history output
      real (kind=dbl_kind), dimension (4,nx_block,ny_block,max_blocks), save :: &
         lont_bounds, & ! longitude of gridbox corners for T point
         latt_bounds, & ! latitude of gridbox corners for T point
         lonu_bounds, & ! longitude of gridbox corners for U point
         latu_bounds    ! latitude of gridbox corners for U point       

      ! geometric quantities used for remapping transport
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         xav  , & ! mean T-cell value of x
         yav  , & ! mean T-cell value of y
         xxav , & ! mean T-cell value of xx
         xyav , & ! mean T-cell value of xy
         yyav , & ! mean T-cell value of yy
         xxxav, & ! mean T-cell value of xxx
         xxyav, & ! mean T-cell value of xxy
         xyyav, & ! mean T-cell value of xyy
         yyyav    ! mean T-cell value of yyy

      real (kind=dbl_kind), &
         dimension (2,2,nx_block,ny_block,max_blocks), save :: &
         mne, & ! matrices used for coordinate transformations in remapping
         mnw, & ! ne = northeast corner, nw = northwest, etc.
         mse, & 
         msw

      ! masks
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         hm     , & ! land/boundary mask, thickness (T-cell)
         bm     , & ! cice block mask (T-cell)
         uvm        ! land/boundary mask, velocity (U-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         ocn_gridcell_frac   ! only relevant for lat-lon grids
                             ! gridcell value of [1 - (land fraction)] (T-cell)

      logical (kind=log_kind), &
         dimension (nx_block,ny_block,max_blocks), save :: &
         tmask  , & ! land/boundary mask, thickness (T-cell)
         umask  , & ! land/boundary mask, velocity (U-cell)
         lmask_n, & ! northern hemisphere mask
         lmask_s    ! southern hemisphere mask

      ! grid dimensions for rectangular grid
      real (kind=dbl_kind), parameter ::  &
         dxrect = 30.e5_dbl_kind   ,&! uniform HTN (cm)
         dyrect = 30.e5_dbl_kind     ! uniform HTE (cm)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), save :: &
         rndex_global       ! global index for local subdomain (dbl)

      real (kind=dbl_kind), private, &
         dimension(nx_block,ny_block,max_blocks) :: &
         work1

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_grid1 - - distribute blocks across processors
!
! !INTERFACE:
!
      subroutine init_grid1 
!
! !DESCRIPTION:
!
! Distribute blocks across processors.  The distribution is optimized
! based on latitude and topography, contained in the ULAT and KMT arrays. 
!
! !REVISION HISTORY:
!
! authors: William Lipscomb and Phil Jones, LANL
!
! !USES:
!
      use ice_broadcast
      use ice_work, only: work_g1, work_g2
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

      integer (kind=int_kind) :: &
         i, j, iblk, &
         fid_grid, &     ! file id for netCDF grid file
         fid_kmt         ! file id for netCDF kmt file

      real (dbl_kind), allocatable, dimension(:,:) :: bStats

      integer(kind=int_kind), allocatable :: WorkPerBlock(:),blockType(:)
      real(kind=dbl_kind), allocatable :: ProbPerBlock(:)

      character (char_len) :: &
         fieldname       ! field name in netCDF file

      !-----------------------------------------------------------------
      ! Get global ULAT and KMT arrays used for block decomposition.
      !-----------------------------------------------------------------

      allocate(work_g1(nx_global,ny_global))
      allocate(work_g2(nx_global,ny_global))

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole'      ) then

         if (trim(grid_format) == 'nc') then

            call ice_open_nc(grid_file,fid_grid)
            call ice_open_nc(kmt_file,fid_kmt)

            fieldname='ulat'
            call ice_read_global_nc(fid_grid,1,fieldname,work_g1,.true.)
            fieldname='kmt'
            call ice_read_global_nc(fid_kmt,1,fieldname,work_g2,.true.)

            if (my_task == master_task) then
               call ice_close_nc(fid_grid)
               call ice_close_nc(fid_kmt)
            endif

         else

            call ice_open(nu_grid,grid_file,64) ! ULAT
            call ice_open(nu_kmt, kmt_file, 32) ! KMT

            call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)  ! ULAT
            call ice_read_global(nu_kmt, 1,work_g2,'ida4',.true.)  ! KMT

            if (my_task == master_task) then
               close (nu_grid)
               close (nu_kmt)
            endif

         endif

      elseif (trim(grid_type) == 'panarctic') then

         call ice_open(nu_grid,grid_file,64) ! ULAT, KMT

         call ice_read_global(nu_grid,1,work_g2,'ida8',.true.)  ! KMT
         call ice_read_global(nu_grid,2,work_g1,'rda8',.true.)  ! ULAT



         if (my_task == master_task) close (nu_grid)

      else   ! rectangular grid

         work_g1(:,:) = 75._dbl_kind/rad_to_deg  ! arbitrary polar latitude
         work_g2(:,:) = c1

      endif

      call broadcast_array(work_g1, master_task)   ! ULAT
      call broadcast_array(work_g2, master_task)   ! KMT

      allocate(WorkPerBlock(nblocks_tot),ProbPerBlock(nblocks_tot),blockType(nblocks_tot))
      allocate(bStats(numCoeff,nblocks_tot))
      call CalcWorkPerBlock(distribution_wght, work_g2,work_g1, &
                WorkPerBlock,ProbPerBlock,blockType,bStats)

!      call abort_ice('init_grid1: after call to CalcWorkPerBlock')
!      stop 'init_grid1: after call to CalcWorkPerBlock'

      call init_domain_distribution(work_g2, work_g1,  &
                WorkPerBlock,ProbPerBlock,blockType,bStats,maxDil)  ! KMT, ULAT
!DBG      print *,'init_grid1: after call to init_domain_distribution'

      deallocate(bStats)
      deallocate(WorkPerBlock,ProbPerBlock,blockType)

      deallocate(work_g1)
      deallocate(work_g2)

      !-----------------------------------------------------------------
      ! write additional domain information
      !-----------------------------------------------------------------

      if (my_task == master_task) then
        write(nu_diag,'(a26,i6)') '  Block size:  nx_block = ',nx_block
        write(nu_diag,'(a26,i6)') '               ny_block = ',ny_block
      endif

      end subroutine init_grid1

!=======================================================================
!BOP
!
! !IROUTINE: init_grid2 - horizontal grid initialization
!
! !INTERFACE:
!
      subroutine init_grid2
!
! !DESCRIPTION:
!
! Horizontal grid initialization:
!
!     U{LAT,LONG} = true {latitude,longitude} of U points
!     HT{N,E} = cell widths on {N,E} sides of T cell
!     ANGLE = angle between local x direction and true east
!     hm = land mask (c1 for ocean points, c0 for land points)
!     D{X,Y}{T,U} = {x,y} spacing centered at {T,U} points
!     T-grid and ghost cell values
!     Various grid quantities needed for dynamics and transport
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_work, only: work_g1
      use ice_exit
      use ice_blocks, only: get_block, block
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
         angle_0, angle_w, angle_s, angle_sw

      logical (kind=log_kind), dimension(nx_block,ny_block,max_blocks):: &
         out_of_range

      type (block) :: &
         this_block           ! block information for current block
      
      !-----------------------------------------------------------------
      ! lat, lon, cell widths, angle, land mask
      !-----------------------------------------------------------------

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole'      ) then
         if (trim(grid_format) == 'nc') then
            call popgrid_nc     ! read POP grid lengths from nc file
         else
            call popgrid        ! read POP grid lengths directly
         endif 
      elseif (trim(grid_type) == 'panarctic') then
         call panarctic_grid    ! pan-Arctic grid
#ifdef CCSMCOUPLED
      elseif (trim(grid_type) == 'latlon') then
         call latlongrid        ! lat lon grid for sequential CCSM (CAM mode)
         return
#endif
      else
         call rectgrid          ! regular rectangular grid
      endif

      !-----------------------------------------------------------------
      ! T-grid cell and U-grid cell quantities
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            tarea(i,j,iblk) = dxt(i,j,iblk)*dyt(i,j,iblk)
            uarea(i,j,iblk) = dxu(i,j,iblk)*dyu(i,j,iblk)
            if (tarea(i,j,iblk) > c0) then
               tarear(i,j,iblk) = c1/tarea(i,j,iblk)
            else
               tarear(i,j,iblk) = c0 ! possible on boundaries
            endif
            if (uarea(i,j,iblk) > c0) then
               uarear(i,j,iblk) = c1/uarea(i,j,iblk)
            else
               uarear(i,j,iblk) = c0 ! possible on boundaries
            endif
            tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

            dxhy(i,j,iblk) = p5*(HTE(i,j,iblk) - HTE(i-1,j,iblk))
            dyhx(i,j,iblk) = p5*(HTN(i,j,iblk) - HTN(i,j-1,iblk))
         enddo
         enddo

         do j = jlo, jhi+1
         do i = ilo, ihi+1
            cyp(i,j,iblk) = (c1p5*HTE(i,j,iblk) - p5*HTE(i-1,j,iblk))
            cxp(i,j,iblk) = (c1p5*HTN(i,j,iblk) - p5*HTN(i,j-1,iblk))
            ! match order of operations in cyp, cxp for tripole grids
            cym(i,j,iblk) = -(c1p5*HTE(i-1,j,iblk) - p5*HTE(i,j,iblk))
            cxm(i,j,iblk) = -(c1p5*HTN(i,j-1,iblk) - p5*HTN(i,j,iblk))
         enddo
         enddo

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Ghost cell updates
      ! On the tripole grid, one must be careful with updates of
      !  quantities that involve a difference of cell lengths.
      ! For example, dyhx and dxhy are cell-centered vector components.
      ! Also note that on the tripole grid, cxp and cxm would swap places,
      !  as would cyp and cym.  These quantities are computed only
      !  in north and east ghost cells (above), not south and west.
      !-----------------------------------------------------------------

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (tarea,              halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (uarea,              halo_info, &
                           field_loc_NEcorner, field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (tarear,             halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (uarear,             halo_info, &
                           field_loc_NEcorner, field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (tinyarea,           halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (dxhy,               halo_info, &
                           field_loc_center,   field_type_vector, &
                           fillValue=c1)
      call ice_HaloUpdate (dyhx,               halo_info, &
                           field_loc_center,   field_type_vector, &
                           fillValue=c1)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! Calculate ANGLET to be compatible with POP ocean model
      ! First, ensure that -pi <= ANGLE <= pi
      !-----------------------------------------------------------------

      out_of_range = .false.
      where (ANGLE < -pi .or. ANGLE > pi) out_of_range = .true.
      if (count(out_of_range) > 0) then
         call abort_ice ('ice: init_grid: ANGLE out of expected range')
      endif

      !-----------------------------------------------------------------
      ! Compute ANGLE on T-grid
      !-----------------------------------------------------------------
      ANGLET = c0

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j,&
      !$OMP                     angle_0,angle_w,angle_s,angle_sw)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            angle_0  = ANGLE(i  ,j  ,iblk) !   w----0
            angle_w  = ANGLE(i-1,j  ,iblk) !   |    |
            angle_s  = ANGLE(i,  j-1,iblk) !   |    |
            angle_sw = ANGLE(i-1,j-1,iblk) !   sw---s

            if ( angle_0 < c0 ) then
               if ( abs(angle_w - angle_0) > pi) &
                        angle_w = angle_w  - pi2
               if ( abs(angle_s - angle_0) > pi) &
                        angle_s = angle_s  - pi2
               if ( abs(angle_sw - angle_0) > pi) &
                        angle_sw = angle_sw - pi2
            endif

            ANGLET(i,j,iblk) = angle_0 * p25 + angle_w * p25 &
                             + angle_s * p25 + angle_sw* p25
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      
      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (ANGLET,           halo_info, &
                           field_loc_center, field_type_angle, &
                           fillValue=c1)
      call ice_timer_stop(timer_bound)

      call makemask          ! velocity mask, hemisphere masks

      call Tlatlon           ! get lat, lon on the T grid

      !----------------------------------------------------------------
      ! Corner coordinates for CF compliant history files
      !----------------------------------------------------------------

      call gridbox_corners

      !-----------------------------------------------------------------
      ! Compute global index (used for unpacking messages from coupler)
      !-----------------------------------------------------------------

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         do j=1,ny_global
         do i=1,nx_global
            work_g1(i,j) = real((j-1)*nx_global + i,kind=dbl_kind)
         enddo
         enddo
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call scatter_global(rndex_global, work_g1,  &
                          master_task,  distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g1)

      end subroutine init_grid2

!=======================================================================
!BOP
!
! !IROUTINE: popgrid - read and set POP displaced pole (or tripole)
!                      grid and land mask
!
! !INTERFACE:
!
      subroutine popgrid
!
! !DESCRIPTION:
!
! POP displaced pole grid and land mask. 
! Grid record number, field and units are: \\
! (1) ULAT  (radians)    \\
! (2) ULON  (radians)    \\
! (3) HTN   (cm)         \\
! (4) HTE   (cm)         \\
! (5) HUS   (cm)         \\
! (6) HUW   (cm)         \\
! (7) ANGLE (radians)   
!
! Land mask record number and field is (1) KMT.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_work, only: work_g1
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      type (block) :: &
         this_block           ! block information for current block

      call ice_open(nu_grid,grid_file,64)
      call ice_open(nu_kmt,kmt_file,32)

      diag = .true.       ! write diagnostic info

      !-----------------------------------------------------------------
      ! topography
      !-----------------------------------------------------------------

      call ice_read(nu_kmt,1,work1,'ida4',diag, &
                    field_loc=field_loc_center, & 
                    field_type=field_type_scalar)

      hm(:,:,:) = c0
      bm(:,:,:) = c0
      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            hm(i,j,iblk) = work1(i,j,iblk)
            if (hm(i,j,iblk) >= c1) hm(i,j,iblk) = c1
            bm(i,j,iblk) = my_task + iblk/1000.0_dbl_kind
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! lat, lon, angle
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1))
      endif

      call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)   ! ULAT
      call gridbox_verts(work_g1,latt_bounds)       
      call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULAT, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      call ice_read_global(nu_grid,2,work_g1,'rda8',.true.)   ! ULON
      call gridbox_verts(work_g1,lont_bounds)       
      call scatter_global(ULON, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULON, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      call ice_read_global(nu_grid,7,work_g1,'rda8',.true.)   ! ANGLE
      call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_angle)

      !-----------------------------------------------------------------
      ! cell dimensions
      ! calculate derived quantities from global arrays to preserve 
      ! information on boundaries
      !-----------------------------------------------------------------

      call ice_read_global(nu_grid,3,work_g1,'rda8',.true.)   ! HTN
      call primary_grid_lengths_HTN(work_g1)                  ! dxu, dxt

      call ice_read_global(nu_grid,4,work_g1,'rda8',.true.)   ! HTE
      call primary_grid_lengths_HTE(work_g1)                  ! dyu, dyt

      deallocate(work_g1)

      if (my_task == master_task) then
         close (nu_grid)
         close (nu_kmt)
      endif

      end subroutine popgrid

!=======================================================================
!BOP
!
! !IROUTINE: popgrid_nc - read and set POP tripole
!                      grid and land mask from netCDF file 
!
! !INTERFACE:
!
      subroutine popgrid_nc

#ifdef ncdf
!
! !DESCRIPTION:
!
! POP displaced pole grid and land mask.
! Grid record number, field and units are: \\
! (1) ULAT  (radians)    \\
! (2) ULON  (radians)    \\
! (3) HTN   (cm)         \\
! (4) HTE   (cm)         \\
! (5) HUS   (cm)         \\
! (6) HUW   (cm)         \\
! (7) ANGLE (radians)
!
! Land mask record number and field is (1) KMT.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
! Revised for netcdf input: Ann Keen, Met Office, May 2007
!
! !USES:
!
      use ice_work, only: work_g1
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi, &     ! beginning and end of physical domain
	 fid_grid, &		! file id for netCDF grid file
	 fid_kmt		! file id for netCDF kmt file

      logical (kind=log_kind) :: diag

      character (char_len) :: &
         fieldname		! field name in netCDF file

      type (block) :: &
         this_block           ! block information for current block

      call ice_open_nc(grid_file,fid_grid)
      call ice_open_nc(kmt_file,fid_kmt)

      diag = .true.       ! write diagnostic info

      !-----------------------------------------------------------------
      ! topography
      !-----------------------------------------------------------------

      fieldname='kmt'
      call ice_read_nc(fid_kmt,1,fieldname,work1,diag, &
                       field_loc=field_loc_center, & 
                       field_type=field_type_scalar)

      hm(:,:,:) = c0
      bm(:,:,:) = c0
      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            hm(i,j,iblk) = work1(i,j,iblk)
            if (hm(i,j,iblk) >= c1) hm(i,j,iblk) = c1
            bm(i,j,iblk) = my_task + iblk/1000.0_dbl_kind
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! lat, lon, angle
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1))
      endif

      fieldname='ulat'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ULAT
      call gridbox_verts(work_g1,latt_bounds)       
      call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULAT, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      fieldname='ulon'
      call ice_read_global_nc(fid_grid,2,fieldname,work_g1,diag) ! ULON
      call gridbox_verts(work_g1,lont_bounds)       
      call scatter_global(ULON, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULON, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      fieldname='angle'
      call ice_read_global_nc(fid_grid,7,fieldname,work_g1,diag) ! ANGLE    
      call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_angle)

      ! fix ANGLE: roundoff error due to single precision
      where (ANGLE >  pi) ANGLE =  pi
      where (ANGLE < -pi) ANGLE = -pi

      !-----------------------------------------------------------------
      ! cell dimensions
      ! calculate derived quantities from global arrays to preserve 
      ! information on boundaries
      !-----------------------------------------------------------------

      fieldname='htn'
      call ice_read_global_nc(fid_grid,3,fieldname,work_g1,diag) ! HTN
      call primary_grid_lengths_HTN(work_g1)                  ! dxu, dxt

      fieldname='hte'
      call ice_read_global_nc(fid_grid,4,fieldname,work_g1,diag) ! HTE
      call primary_grid_lengths_HTE(work_g1)                  ! dyu, dyt

      deallocate(work_g1)

      if (my_task == master_task) then
         call ice_close_nc(fid_grid)
         call ice_close_nc(fid_kmt)
      endif

#endif
      end subroutine popgrid_nc

!=======================================================================
!BOP
!
! !IROUTINE: panarctic_grid - read and set Pan-Arctic grid and land mask
!
! !INTERFACE:
!
      subroutine panarctic_grid
!
! !DESCRIPTION:
!
! Pan-Arctic grid and mask developed by Wieslaw Maslowski
!
! !REVISION HISTORY:
!
! authors: Wieslaw Maslowki, Naval Postgraduate School (based on popgrid)
!          William H. Lipscomb, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_work, only: work_g1
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      !-----------------------------------------------------------------
      ! 
      ! PIPS rotated spherical grid and land mask
      !      rec no.         field         units
      !      -------         -----         -----
      !   land mask
      !         1             KMT         
      !   grid
      !         2            ULAT         radians
      !         3            ULON         radians
      !         4             HTN           cm
      !         5             HTE           cm
      !         6             HUS           cm
      !         7             HUW           cm
      !         8            ANGLE        radians
      !
      ! NOTE: There is no separate kmt file.  Land mask is part of grid file.
      !-----------------------------------------------------------------

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      type (block) :: &
         this_block           ! block information for current block

      call ice_open(nu_grid,grid_file,64)

      diag = .true.       ! write diagnostic info

      if (my_task == master_task) &
           write (nu_diag,*) '** Reading pan-Arctic grid **'

      !-----------------------------------------------------------------
      ! topography
      !-----------------------------------------------------------------

      call ice_read(nu_grid,1,work1,'ida8',diag, &
                    field_loc=field_loc_center, & 
                    field_type=field_type_scalar)

      hm(:,:,:) = c0
      bm(:,:,:) = c0
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            hm(i,j,iblk) = work1(i,j,iblk)
            if (hm(i,j,iblk) >= c1) hm(i,j,iblk) = c1
            bm(i,j,iblk) = my_task + iblk/1000.0_dbl_kind
         enddo
         enddo
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! lat, lon, angle
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1))
      endif

      call ice_read_global(nu_grid,2,work_g1,'rda8',.true.)   ! ULAT
      call gridbox_verts(work_g1,latt_bounds)       
      call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULAT, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      call ice_read_global(nu_grid,3,work_g1,'rda8',.true.)   ! ULON
      call gridbox_verts(work_g1,lont_bounds)       
      call scatter_global(ULON, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULON, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      call ice_read_global(nu_grid,8,work_g1,'rda8',.true.)   ! ANGLE
      call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_angle)

      !-----------------------------------------------------------------
      ! cell dimensions
      ! calculate derived quantities from global arrays to preserve 
      ! information on boundaries
      !-----------------------------------------------------------------

      call ice_read_global(nu_grid,4,work_g1,'rda8',.true.)   ! HTN
      call primary_grid_lengths_HTN(work_g1)                  ! dxu, dxt

      call ice_read_global(nu_grid,5,work_g1,'rda8',.true.)   ! HTE
      call primary_grid_lengths_HTE(work_g1)                  ! dyu, dyt

      deallocate(work_g1)

      if (my_task == master_task) close (nu_grid)

      end subroutine panarctic_grid
#ifdef CCSMCOUPLED
!=======================================================================
!BOP
!
! !IROUTINE: latlongrid- lat and lon grid for coupling to standalone CAM
!
! !INTERFACE:
!
      subroutine latlongrid
!
! !DESCRIPTION:
!
! Read in kmt file that matches CAM lat-lon grid and has single column 
! functionality
!
! !REVISION HISTORY:
!
! author: Mariana Vertenstein
! 2007: Elizabeth Hunke upgraded to netcdf90 and cice ncdf calls
!
! !USES:
!
!     use ice_boundary
      use ice_domain_size
      use ice_scam, only : scmlat, scmlon, single_column
      use netcdf
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk    
      
      integer (kind=int_kind) :: &
         ni, nj, ncid, dimid, varid, ier

      character (len=char_len) :: &
         subname='latlongrid' ! subroutine name

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
           closelat, &        ! Single-column latitude value
           closelon, &        ! Single-column longitude value
           closelatidx, &     ! Single-column latitude index to retrieve
           closelonidx        ! Single-column longitude index to retrieve

      integer (kind=int_kind) :: &
           start(2), &        ! Start index to read in
           count(2)           ! Number of points to read in

      integer (kind=int_kind) :: &
           start3(3), &        ! Start index to read in
           count3(3)           ! Number of points to read in

      integer (kind=int_kind) :: &
        status                ! status flag

      real (kind=dbl_kind), allocatable :: &
           lats(:),lons(:),pos_lons(:), glob_grid(:,:)  ! temporaries 

      real (kind=dbl_kind) :: &
         pos_scmlon,&         ! temporary
         scamdata,&           ! temporary
         minpoint,&
         testpoint

      logical :: se_grid      ! determine whether file format is for 
                              !   SE dynamical core grid or not

      !-----------------------------------------------------------------
      ! - kmt file is actually clm fractional land file
      ! - Determine consistency of dimensions
      ! - Read in lon/lat centers in degrees from kmt file
      ! - Read in ocean from "kmt" file (1 for ocean, 0 for land)
      !-----------------------------------------------------------------

      ! Determine dimension of domain file and check for consistency

      if (my_task == master_task) then
         call ice_open_nc(kmt_file, ncid)

         status = nf90_inq_dimid (ncid, 'ni', dimid)
         status = nf90_inquire_dimension(ncid, dimid, len=ni)
         status = nf90_inq_dimid (ncid, 'nj', dimid)
         status = nf90_inquire_dimension(ncid, dimid, len=nj)
      end if

      ! Determine start/count to read in for either single column or global lat-lon grid
      ! If single_column, then assume that only master_task is used since there is only one task

      if (single_column) then
         ! Check for consistency
         if (my_task == master_task) then
            if ((nx_global /= 1).or. (ny_global /= 1)) then
               write(nu_diag,*) 'Because you have selected the column model flag'
               write(nu_diag,*) 'Please set nx_global=ny_global=1 in file'
               write(nu_diag,*) 'ice_domain_size.F and recompile'
               call abort_ice ('latlongrid: check nx_global, ny_global')
            endif
         end if

         ! Are we dealing with a file for SE grid?
         se_grid=.false.
         if (nj .eq. 1 .and. ni .gt. 1) then
           se_grid=.true.
         endif

         ! Read in domain file for single column
         if (se_grid) then
           allocate(lats(ni))
         else
           allocate(lats(nj))
         endif

         allocate(lons(ni))
         allocate(pos_lons(ni))
         allocate(glob_grid(ni,nj))

         start3=(/1,1,1/)
         count3=(/ni,nj,1/)
         call check_ret(nf90_inq_varid(ncid, 'xc' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, glob_grid, start3, count3), subname)
         do i = 1,ni
            lons(i) = glob_grid(i,1)
         end do

         call check_ret(nf90_inq_varid(ncid, 'yc' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, glob_grid, start3, count3), subname)

         if (se_grid) then
           do i=1,ni
             lats(i) = glob_grid(i,1)
           end do
         else
           do j=1,nj
             lats(j) = glob_grid(1,j)
           enddo
         endif

         ! convert lons array and scmlon to 0,360 and find index of value closest to 0
         ! and obtain single-column longitude/latitude indices to retrieve
         
         pos_lons(:)= mod(lons(:) + 360._dbl_kind,360._dbl_kind)
         pos_scmlon = mod(scmlon  + 360._dbl_kind,360._dbl_kind)

         if (se_grid) then

           minpoint=1000.0
           do i=1,ni
             testpoint=abs(pos_lons(i)-pos_scmlon)+abs(lats(i)-scmlat)
             if (testpoint .lt. minpoint) then
               minpoint=testpoint
               start(1)=i
               start(2)=1
             endif
           enddo

         else

           start(1) = (MINLOC(abs(pos_lons-pos_scmlon),dim=1))
           start(2) = (MINLOC(abs(lats    -scmlat    ),dim=1))

         endif

         deallocate(lats)
         deallocate(lons)
         deallocate(pos_lons)
         deallocate(glob_grid)

         call check_ret(nf90_inq_varid(ncid, 'xc' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, scamdata, start), subname)
         TLON = scamdata
         call check_ret(nf90_inq_varid(ncid, 'yc' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, scamdata, start), subname)
         TLAT = scamdata
         call check_ret(nf90_inq_varid(ncid, 'area' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, scamdata, start), subname)
         tarea = scamdata
         call check_ret(nf90_inq_varid(ncid, 'mask' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, scamdata, start), subname)
         hm = scamdata
         call check_ret(nf90_inq_varid(ncid, 'frac' , varid), subname)
         call check_ret(nf90_get_var(ncid, varid, scamdata, start), subname)
         ocn_gridcell_frac = scamdata
      else
         ! Check for consistency
         if (my_task == master_task) then
            if (nx_global /= ni .and. ny_global /= nj) then
              call abort_ice ('latlongrid: ni,nj not equal to nx_global,ny_global')
            end if
         end if

         ! Read in domain file for global lat-lon grid
         call ice_read_nc(ncid, 1, 'xc'  , TLON             , diag=.true.)
         call ice_read_nc(ncid, 1, 'yc'  , TLAT             , diag=.true.)
         call ice_read_nc(ncid, 1, 'area', tarea            , diag=.true., &
            field_loc=field_loc_center,field_type=field_type_scalar)
         call ice_read_nc(ncid, 1, 'mask', hm               , diag=.true.)
         call ice_read_nc(ncid, 1, 'frac', ocn_gridcell_frac, diag=.true.)
      end if

      if (my_task == master_task) then
         call ice_close_nc(ncid)
      end if

     !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1,nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            ! Convert from degrees to radians
            TLON(i,j,iblk) = pi*TLON(i,j,iblk)/180._dbl_kind

            ! Convert from degrees to radians
            TLAT(i,j,iblk) = pi*TLAT(i,j,iblk)/180._dbl_kind

            ! Convert from radians^2 to m^2
            ! (area in domain file is in radians^2 and tarea is in m^2)
            tarea(i,j,iblk) = tarea(i,j,iblk) * (radius*radius)
            bm(i,j,iblk) = my_task + iblk/1000.0_dbl_kind
         end do
         end do
      end do
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      ! The U grid (velocity) is not used when run with sequential CAM
      ! because we only use thermodynamic sea ice.  However, ULAT is used
      ! in the default initialization of CICE so we calculate it here as 
      ! a "dummy" so that CICE will initialize with ice.  If a no ice
      ! initialization is OK (or desired) this can be commented out and
      ! ULAT will remain 0 as specified above.  ULAT is located at the
      ! NE corner of the grid cell, TLAT at the center, so here ULAT is
      ! hacked by adding half the latitudinal spacing (in radians) to
      ! TLAT.
      !-----------------------------------------------------------------

     !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            if (ny_global == 1) then
               uarea(i,j,iblk)  = tarea(i,j,  iblk)
            else
               uarea(i,j,iblk)  = p25*  &
                                 (tarea(i,j,  iblk) + tarea(i+1,j,  iblk) &
                                + tarea(i,j+1,iblk) + tarea(i+1,j+1,iblk))
            endif
            tarear(i,j,iblk)   = c1/tarea(i,j,iblk)
            uarear(i,j,iblk)   = c1/uarea(i,j,iblk)
            tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

            if (single_column) then
               ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/nj)  
            else
               if (ny_global == 1) then
                  ULAT  (i,j,iblk) = TLAT(i,j,iblk)
               else
                  ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/ny_global)  
               endif
            endif
            ULON  (i,j,iblk) = c0
            ANGLE (i,j,iblk) = c0                             

            ANGLET(i,j,iblk) = c0                             
            HTN   (i,j,iblk) = 1.e36_dbl_kind
            HTE   (i,j,iblk) = 1.e36_dbl_kind
            dxt   (i,j,iblk) = 1.e36_dbl_kind
            dyt   (i,j,iblk) = 1.e36_dbl_kind
            dxu   (i,j,iblk) = 1.e36_dbl_kind
            dyu   (i,j,iblk) = 1.e36_dbl_kind
            dxhy  (i,j,iblk) = 1.e36_dbl_kind
            dyhx  (i,j,iblk) = 1.e36_dbl_kind
            cyp   (i,j,iblk) = 1.e36_dbl_kind
            cxp   (i,j,iblk) = 1.e36_dbl_kind
            cym   (i,j,iblk) = 1.e36_dbl_kind
            cxm   (i,j,iblk) = 1.e36_dbl_kind
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call makemask

      end subroutine latlongrid

!=======================================================================
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
      subroutine check_ret(ret, calling)
!
! !DESCRIPTION:
!     Check return status from netcdf call
!
! !USES:
        use netcdf
! !ARGUMENTS:
        implicit none
        integer, intent(in) :: ret
        character(len=*) :: calling
!
! !REVISION HISTORY:
! author: Mariana Vertenstein
! 2007: Elizabeth Hunke upgraded to netcdf90
!
!EOP
!
      if (ret /= NF90_NOERR) then
         write(nu_diag,*)'netcdf error from ',trim(calling)
         write(nu_diag,*)'netcdf strerror = ',trim(NF90_STRERROR(ret))
         call abort_ice('ice ice_grid: netcdf check_ret error')
      end if
        
      end subroutine check_ret

#endif      
!=======================================================================
!BOP
!
! !IROUTINE: rectgrid - regular rectangular grid and mask
!
! !INTERFACE:
!
      subroutine rectgrid
!
! !DESCRIPTION:
!
! Regular rectangular grid and mask
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_work, only: work_g1
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         imid, jmid

      real (kind=dbl_kind) :: length

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            ANGLE(i,j,iblk) = c0              ! "square with the world"
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1))
      endif

      ! Weddell Sea
      ! lower left corner of grid is 55W, 75S

      if (my_task == master_task) then
         work_g1 = c0
         length = dxrect*cm_to_m/radius*rad_to_deg
         work_g1(1,:) = -55._dbl_kind
         do j = 1, ny_global
         do i = 2, nx_global
            work_g1(i,j) = work_g1(i-1,j) + length   ! ULON
         enddo
         enddo
         work_g1(:,:) = work_g1(:,:) / rad_to_deg
      endif
      call scatter_global(ULON, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      if (my_task == master_task) then
         work_g1 = c0
         length = dyrect*cm_to_m/radius*rad_to_deg
         work_g1(:,1) = -75._dbl_kind
         do i = 1, nx_global
         do j = 2, ny_global
            work_g1(i,j) = work_g1(i,j-1) + length   ! ULAT
         enddo
         enddo
         work_g1(:,:) = work_g1(:,:) / rad_to_deg
      endif
      call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 1, nx_global
            work_g1(i,j) = dxrect             ! HTN
         enddo
         enddo
      endif
      call primary_grid_lengths_HTN(work_g1)  ! dxu, dxt

      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 1, nx_global
            work_g1(i,j) = dyrect             ! HTE
         enddo
         enddo
      endif
      call primary_grid_lengths_HTE(work_g1)  ! dyu, dyt

      !-----------------------------------------------------------------
      ! Construct T-cell land mask
      ! Keyed on ew_boundary_type; ns_boundary_type should be 'open'.
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         work_g1(:,:) = c0      ! initialize hm as land

         if (trim(ew_boundary_type) == 'cyclic') then

            do j = 3,ny_global-2      ! closed top and bottom
            do i = 1,nx_global        ! open sides
               work_g1(i,j) = c1    ! NOTE nx_global > 5
            enddo
            enddo

         elseif (trim(ew_boundary_type) == 'closed') then

            do j = 3,ny_global-2      ! closed top and bottom
            do i = 3,nx_global-2      ! closed sides
               work_g1(i,j) = c1    ! NOTE nx_global, ny_global > 5
            enddo
            enddo

         elseif (trim(ew_boundary_type) == 'open') then

            ! land in the upper left and lower right corners,
            ! otherwise open boundaries
            imid = aint(real(nx_global)/c2,kind=int_kind)
            jmid = aint(real(ny_global)/c2,kind=int_kind)

            do j = 3,ny_global-2
            do i = 3,nx_global-2
               work_g1(i,j) = c1    ! open central domain
            enddo
            enddo

            do j = 1, jmid+2
            do i = 1, imid+2
               work_g1(i,j) = c1    ! open lower left corner
            enddo
            enddo

            do j = jmid-2, ny_global
            do i = imid-2, nx_global
               work_g1(i,j) = c1    ! open upper right corner
            enddo
            enddo

         endif
      endif

      call scatter_global(hm, work_g1, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g1)

      end subroutine rectgrid

!=======================================================================
!BOP
!
! !IROUTINE: primary_grid_lengths_HTN
!
! !INTERFACE:
!
      subroutine primary_grid_lengths_HTN(work_g)
!
! !DESCRIPTION:
!
! Calculate dxu and dxt from HTN on the global grid, to preserve
! ghost cell and/or land values that might otherwise be lost. Scatter
! dxu, dxt and HTN to all processors.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_work, only: work_g2
!
! !INPUT/OUTPUT PARAMETERS:
!
! work_g is the global array holding HTN.
!
!EOP

      real (kind=dbl_kind), dimension(:,:) :: work_g

      integer (kind=int_kind) :: &
         i, j, &
         ip1     ! i+1

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global))
      else
         allocate(work_g2(1,1))
      endif

      if (my_task == master_task) then
      do j = 1, ny_global
      do i = 1, nx_global
         work_g(i,j) = work_g(i,j) * cm_to_m                ! HTN
      enddo
      enddo
      do j = 1, ny_global
      do i = 1, nx_global
         ! assume cyclic; noncyclic will be handled during scatter
         ip1 = i+1
         if (i == nx_global) ip1 = 1
         work_g2(i,j) = p5*(work_g(i,j) + work_g(ip1,j))    ! dxu
      enddo
      enddo
      endif
      call scatter_global(HTN, work_g, master_task, distrb_info, &
                          field_loc_Nface, field_type_scalar)
      call scatter_global(dxu, work_g2, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      if (my_task == master_task) then
      do j = 2, ny_global
         do i = 1, nx_global
            work_g2(i,j) = p5*(work_g(i,j) + work_g(i,j-1)) ! dxt
         enddo
      enddo
      ! extrapolate to obtain dxt along j=1
      do i = 1, nx_global
         work_g2(i,1) = c2*work_g(i,2) - work_g(i,3) ! dxt
      enddo
      endif
      call scatter_global(dxt, work_g2, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g2)

      end subroutine primary_grid_lengths_HTN

!=======================================================================
!BOP
!
! !IROUTINE: primary_grid_lengths_HTE
!
! !INTERFACE:
!
      subroutine primary_grid_lengths_HTE(work_g)
!
! !DESCRIPTION:
!
! Calculate dyu and dyt from HTE on the global grid, to preserve
! ghost cell and/or land values that might otherwise be lost. Scatter
! dyu, dyt and HTE to all processors.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_work, only: work_g2
!
! !INPUT/OUTPUT PARAMETERS:
!
! work_g is the global array holding HTE.
!
!EOP

      real (kind=dbl_kind), dimension(:,:) :: work_g

      integer (kind=int_kind) :: &
         i, j, nyg, &
         im1     ! i-1

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global))
      else
         allocate(work_g2(1,1))
      endif

      if (my_task == master_task) then
      do j = 1, ny_global
      do i = 1, nx_global
         work_g(i,j) = work_g(i,j) * cm_to_m                ! HTE
      enddo
      enddo
      do j = 1, ny_global-1
         do i = 1, nx_global
            work_g2(i,j) = p5*(work_g(i,j) + work_g(i,j+1)) ! dyu
         enddo
      enddo
      ! extrapolate to obtain dyu along j=ny_global
      ! workaround for intel compiler
      nyg = ny_global
      do i = 1, nx_global
         work_g2(i,nyg) = c2*work_g(i,nyg-1) &
                           - work_g(i,nyg-2) ! dyu
      enddo
      endif
      call scatter_global(HTE, work_g, master_task, distrb_info, &
                          field_loc_Eface, field_type_scalar)
      call scatter_global(dyu, work_g2, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      if (my_task == master_task) then
      do j = 1, ny_global
      do i = 1, nx_global
         ! assume cyclic; noncyclic will be handled during scatter
         im1 = i-1
         if (i == 1) im1 = nx_global 
         work_g2(i,j) = p5*(work_g(i,j) + work_g(im1,j))    ! dyt
      enddo
      enddo
      endif
      call scatter_global(dyt, work_g2, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g2)

      end subroutine primary_grid_lengths_HTE

!=======================================================================
!BOP
!
! !IROUTINE: makemask - makes logical land masks (T,U) and hemispheric masks
!
! !INTERFACE:
!
      subroutine makemask
!
! !DESCRIPTION:
!
! Sets the boundary values for the T cell land mask (hm) and
! makes the logical land masks for T and U cells (tmask, umask).
! Also creates hemisphere masks (mask-n northern, mask-s southern)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (hm,               halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_HaloUpdate (bm,               halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! construct T-cell and U-cell masks
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            if (ny_global == 1) then
               uvm(i,j,iblk) =      hm(i,j,  iblk)
            else
               uvm(i,j,iblk) = min (hm(i,j,  iblk), hm(i+1,j,  iblk), &
                                    hm(i,j+1,iblk), hm(i+1,j+1,iblk))
            endif
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (uvm,                halo_info, &
                           field_loc_NEcorner, field_type_scalar)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            tmask(i,j,iblk) = .false.
            umask(i,j,iblk) = .false.
            if ( hm(i,j,iblk) > p5) tmask(i,j,iblk) = .true.
            if (uvm(i,j,iblk) > p5) umask(i,j,iblk) = .true.
         enddo
         enddo

      !-----------------------------------------------------------------
      ! create hemisphere masks
      !-----------------------------------------------------------------

         lmask_n(:,:,iblk) = .false.
         lmask_s(:,:,iblk) = .false.

         tarean(:,:,iblk) = c0
         tareas(:,:,iblk) = c0

         do j = 1, ny_block
         do i = 1, nx_block

            if (ULAT(i,j,iblk) >= -puny) lmask_n(i,j,iblk) = .true. ! N. Hem.
            if (ULAT(i,j,iblk) <  -puny) lmask_s(i,j,iblk) = .true. ! S. Hem.

            ! N hemisphere area mask (m^2)
            if (lmask_n(i,j,iblk)) tarean(i,j,iblk) = tarea(i,j,iblk) &
                                                    * hm(i,j,iblk)

            ! S hemisphere area mask (m^2)
            if (lmask_s(i,j,iblk)) tareas(i,j,iblk) = tarea(i,j,iblk) &
                                                    * hm(i,j,iblk)

         enddo
         enddo

      enddo  ! iblk
      !$OMP END PARALLEL DO

      end subroutine makemask

!=======================================================================
!BOP
!
! !IROUTINE: Tlatlon - initializes latitude and longitudes on T grid
!
! !INTERFACE:
!
      subroutine Tlatlon
!
! !DESCRIPTION:
!
! Initializes latitude and longitude on T grid
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL; code originally based on POP grid
! generation routine
!
! !USES:
!
      use ice_domain_size
      use ice_global_reductions, only: global_minval, global_maxval
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      save 

      integer (kind=int_kind) :: &
           i, j, iblk       , & ! horizontal indices
           ig, jg           , & ! global horizontal indices
           im1              , & ! ig - 1
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
           z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da

      type (block) :: &
           this_block           ! block information for current block

      TLAT(:,:,:) = c0
      TLON(:,:,:) = c0

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j, &
      !$OMP                    z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4, &
      !$OMP                     tx,ty,tz,da)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            z1 = cos(ULAT(i-1,j-1,iblk))
            x1 = cos(ULON(i-1,j-1,iblk))*z1
            y1 = sin(ULON(i-1,j-1,iblk))*z1
            z1 = sin(ULAT(i-1,j-1,iblk))

            z2 = cos(ULAT(i,j-1,iblk))
            x2 = cos(ULON(i,j-1,iblk))*z2
            y2 = sin(ULON(i,j-1,iblk))*z2
            z2 = sin(ULAT(i,j-1,iblk))

            z3 = cos(ULAT(i-1,j,iblk))
            x3 = cos(ULON(i-1,j,iblk))*z3
            y3 = sin(ULON(i-1,j,iblk))*z3
            z3 = sin(ULAT(i-1,j,iblk))

            z4 = cos(ULAT(i,j,iblk))
            x4 = cos(ULON(i,j,iblk))*z4
            y4 = sin(ULON(i,j,iblk))*z4
            z4 = sin(ULAT(i,j,iblk))

            tx = (x1+x2+x3+x4)/c4
            ty = (y1+y2+y3+y4)/c4
            tz = (z1+z2+z3+z4)/c4
            da = sqrt(tx**2+ty**2+tz**2)

            tz = tz/da

            ! TLON in radians East
            TLON(i,j,iblk) = c0
            if (tx /= c0 .or. ty /= c0) TLON(i,j,iblk) = atan2(ty,tx)

            ! TLAT in radians North
            TLAT(i,j,iblk) = asin(tz)
            
         enddo                  ! i
         enddo                  ! j         
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (TLON,             halo_info, &
                           field_loc_center, field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (TLAT,             halo_info, &
                           field_loc_center, field_type_scalar, &
                           fillValue=c1)
      call ice_HaloExtrapolate(TLON, distrb_info, &
                               ew_boundary_type, ns_boundary_type)
      call ice_HaloExtrapolate(TLAT, distrb_info, &
                               ew_boundary_type, ns_boundary_type)
      call ice_timer_stop(timer_bound)

      x1 = global_minval(TLON, distrb_info, tmask)
      x2 = global_maxval(TLON, distrb_info, tmask)
      x3 = global_minval(TLAT, distrb_info, tmask)
      x4 = global_maxval(TLAT, distrb_info, tmask)

      y1 = global_minval(ULON, distrb_info, umask)
      y2 = global_maxval(ULON, distrb_info, umask)
      y3 = global_minval(ULAT, distrb_info, umask)
      y4 = global_maxval(ULAT, distrb_info, umask)

      if (my_task==master_task) then
         write(nu_diag,*) ' '
         write(nu_diag,*) 'min/max ULON:', y1*rad_to_deg, y2*rad_to_deg
         write(nu_diag,*) 'min/max TLON:', x1*rad_to_deg, x2*rad_to_deg
         write(nu_diag,*) 'min/max ULAT:', y3*rad_to_deg, y4*rad_to_deg
         write(nu_diag,*) 'min/max TLAT:', x3*rad_to_deg, x4*rad_to_deg
      endif                     ! my_task

      end subroutine Tlatlon

!=======================================================================
!BOP
!
! !IROUTINE: t2ugrid_vector - transfer vector from T-cells to U-cells
!
! !INTERFACE:
!
      subroutine t2ugrid_vector (work)
!
! !DESCRIPTION:
!
! Transfer vector component from T-cell centers to U-cell centers.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
           intent(inout) :: & 
           work
      
      integer (int_kind) :: iblk
!
!EOP
!
      work1(:,:,:) = work(:,:,:)
!      do iblk = 1,nblocks
!         work1(:,:,iblk) = work(:,:,iblk)
!      enddo

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (work1,            halo_info, &
                           field_loc_center, field_type_vector)
      call ice_timer_stop(timer_bound)

      call to_ugrid(work1,work)

      end subroutine t2ugrid_vector

!=======================================================================
!BOP
!
! !IROUTINE: to_ugrid - shift from T-cell to U-cell midpoints
!
! !INTERFACE:
!
      subroutine to_ugrid(work1,work2)
!
! !DESCRIPTION:
!
! Shifts quantities from the T-cell midpoint (work1) to the U-cell
! midpoint (work2)
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
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
         work1(nx_block,ny_block,max_blocks)

      real (kind=dbl_kind), intent(out) :: &
         work2(nx_block,ny_block,max_blocks)

      type (block) :: &
         this_block           ! block information for current block
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      work2(:,:,:) = c0

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            work2(i,j,iblk) = p25 * &
                              (work1(i,  j,  iblk)*tarea(i,  j,  iblk)  &
                             + work1(i+1,j,  iblk)*tarea(i+1,j,  iblk)  &
                             + work1(i,  j+1,iblk)*tarea(i,  j+1,iblk)  &
                             + work1(i+1,j+1,iblk)*tarea(i+1,j+1,iblk)) &
                             / uarea(i,  j,  iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine to_ugrid

!=======================================================================
!BOP
!
! !IROUTINE: u2tgrid_vector - transfer vector from U-cells to T-cells
!
! !INTERFACE:
!
      subroutine u2tgrid_vector (work)
!
! !DESCRIPTION:
!
! Transfer from U-cell centers to T-cell centers. Writes work into
! another array that has ghost cells
! NOTE: Input array is dimensioned only over physical cells.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), &
         intent(inout) :: &
         work
!
!EOP
!
      work1(:,:,:) = work(:,:,:)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (work1,            halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call ice_timer_stop(timer_bound)

      call to_tgrid(work1,work)

      end subroutine u2tgrid_vector

!=======================================================================
!BOP
!
! !IROUTINE: to_tgrid - shifts array from U-cell to T-cell midpoints
!
! !INTERFACE:
!
      subroutine to_tgrid(work1, work2)
!
! !DESCRIPTION:
!
! Shifts quantities from the U-cell midpoint (work1) to the T-cell
! midpoint (work2)
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind) :: work1(nx_block,ny_block,max_blocks), &
                              work2(nx_block,ny_block,max_blocks)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block
      
      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            work2(i,j,iblk) = p25 *  &
                             (work1(i,  j  ,iblk) * uarea(i,  j,  iblk)  &
                            + work1(i-1,j  ,iblk) * uarea(i-1,j,  iblk)  &
                            + work1(i,  j-1,iblk) * uarea(i,  j-1,iblk)  & 
                            + work1(i-1,j-1,iblk) * uarea(i-1,j-1,iblk)) &
                            / tarea(i,  j,  iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine to_tgrid

!=======================================================================
!BOP
!
! !IROUTINE: sinlat - calculates sin of latitudes
!
! !INTERFACE:
!
      subroutine sinlat(a, k)
!
! !DESCRIPTION:
!
! Calculates the sin of latitudes on the grid based on ny_global and
! using code from CAM (/control/gauaw_mod.F90)  In CAM, gauaw_mod.F90 is
! only used to calculate the sin of latitudes (and thence latitude) if
! the dynamical core is set to eul.  If using one of the other dynamical
! cores and coupling to stand alone CAM the latitudes calculated in this
! way may not match the grid from CAM.
!
! !REVISION HISTORY:
!
! author: Jacob Sewall
!
! !USES:
      use ice_exit

!
! !INPUT/OUTPUT PARAMETERS:
!

      real (kind=dbl_kind), dimension(ny_global), intent(out) :: & 
           a            ! sin of latitudes

      integer (kind=int_kind), intent(in) :: &
           k            ! number of latitudes (ny_global)
!
!EOP
!
      real (kind=dbl_kind), dimension(k) :: &
           sinlats      ! sine of latitudes

      real (kind=dbl_kind) :: &
           eps      , & ! convergence criterion
           c        , & ! constant combination
           fk       , & ! real k
           xz       , & ! abscissa estimate
           pkm1     , & ! |
           pkm2     , & ! |-polynomials
           pkmrk    , & ! |
           pk       , & ! |
           sp       , & ! current iteration latitude increment
           avsp     , & ! |sp|
           fn       , & ! real n
           avsp_prev, &
           testdiff

      real (kind=dbl_kind), parameter :: &
           eps27 = 1.e-27_dbl_kind

      integer (kind=int_kind) :: &
           n, l     , &
           iter     , & ! iteration counter
           is       , & ! latitude index
           kk           ! k/2 (number of latitudes in a hemisphere)

!
!----------------------------------------------------------------------
!

      c = (c1-(c2/pi)**2)*p25
      fk = k
      kk = k/2
      call bsslzr(sinlats,kk)
      do is=1,kk
        xz = cos(sinlats(is)/sqrt((fk+p5)**2+c))
!
! This is the first approximation to xz
!
        iter = 0
        avsp = c100     !initialize avsp to a very large number
10      continue
        avsp_prev = avsp
        pkm2 = c1
        pkm1 = xz
        iter = iter + 1
        if (iter > 100) then  ! Error exit
           call abort_ice('SINLAT: no convergence in 100 iterations')
        end if
!
! Computation of the legendre polynomial
!
        do n=2,k
          fn = n
          pk = ((c2*fn-1._dbl_kind)*xz*pkm1-(fn-c1)*pkm2)/fn
          pkm2 = pkm1
          pkm1 = pk
        enddo
        pkm1 = pkm2
        pkmrk = (fk*(pkm1-xz*pk))/(c1-xz**2)
        sp = pk/pkmrk
        xz = xz - sp
        avsp = abs(sp)
        testdiff = avsp_prev - avsp
        if (testdiff > eps27) go to 10
        sinlats(is) = xz
      end do
!
! Complete the sets of abscissas and weights, using the symmetry.
! Also note truncation from real(r8) to real*8
!
      do n=1,kk
        l = k + 1 - n
        a(n) = sinlats(n)
        a(l) = -sinlats(n)
      end do
 
      end subroutine sinlat

!=======================================================================
!BOP
!
! !IROUTINE: bsslzr - Return n zeros (if n<50) of the Bessel function
!
! !INTERFACE:
!
      subroutine bsslzr(bes, n)
!
! !DESCRIPTION:
!
!
! Return n zeros (or if n>50, approximate zeros), of the Bessel function
! j0,in the array bes. The first 50 zeros will be given exactly, and the
! remaining zeros are computed by extrapolation,and therefore not exact.
!
! Modified 1/23/97 by Jim Rosinski to use real*16 arithmetic
! placed in ice_grid.F 4/26/2006 by Jacob Sewall

! !REVISION HISTORY:
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, D. Williamson, August 1992
! Reviewed:          J. Hack, D. Williamson, April 1996
!
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!

      integer (kind=int_kind), intent(in) :: &
           n          ! number of latitudes in hemisphere (ny_global/2)

      real (kind=dbl_kind), dimension(n), intent(inout) :: & 
           bes        ! sin of latitudes
!
!EOP
!
!----------------------------------------
! Local Variables
!----------------------------------------

      integer (kind=int_kind) :: &
         nn, j
     
      real (kind=dbl_kind), dimension(50) :: bz
      
      save bz
!
!-----------------------------------------
! Local Workspace
!-----------------------------------------

      data bz/ &
        2.4048255577_dbl_kind,   5.5200781103_dbl_kind,   8.6537279129_dbl_kind, &
       11.7915344391_dbl_kind,  14.9309177086_dbl_kind,  18.0710639679_dbl_kind, &
       21.2116366299_dbl_kind,  24.3524715308_dbl_kind,  27.4934791320_dbl_kind, &
       30.6346064684_dbl_kind,  33.7758202136_dbl_kind,  36.9170983537_dbl_kind, &
       40.0584257646_dbl_kind,  43.1997917132_dbl_kind,  46.3411883717_dbl_kind, &
       49.4826098974_dbl_kind,  52.6240518411_dbl_kind,  55.7655107550_dbl_kind, &
       58.9069839261_dbl_kind,  62.0484691902_dbl_kind,  65.1899648002_dbl_kind, &
       68.3314693299_dbl_kind,  71.4729816036_dbl_kind,  74.6145006437_dbl_kind, &
       77.7560256304_dbl_kind,  80.8975558711_dbl_kind,  84.0390907769_dbl_kind, &
       87.1806298436_dbl_kind,  90.3221726372_dbl_kind,  93.4637187819_dbl_kind, &
       96.6052679510_dbl_kind,  99.7468198587_dbl_kind, 102.8883742542_dbl_kind, &
      106.0299309165_dbl_kind, 109.1714896498_dbl_kind, 112.3130502805_dbl_kind, &
      115.4546126537_dbl_kind, 118.5961766309_dbl_kind, 121.7377420880_dbl_kind, &
      124.8793089132_dbl_kind, 128.0208770059_dbl_kind, 131.1624462752_dbl_kind, &
      134.3040166383_dbl_kind, 137.4455880203_dbl_kind, 140.5871603528_dbl_kind, &
      143.7287335737_dbl_kind, 146.8703076258_dbl_kind, 150.0118824570_dbl_kind, &
      153.1534580192_dbl_kind, 156.2950342685_dbl_kind/  

      nn = n 
      if (n > 50) then 
         bes(50) = bz(50) 
         do j = 51, n 
            bes(j) = bes(j-1) + pi 
         end do 
         nn = 49 
      endif 
      bes(:nn) = bz(:nn) 

      end subroutine bsslzr 

!=======================================================================
! The following code is used for obtaining the coordinates of the grid
! vertices for CF-compliant netCDF history output. Approximate!
!=======================================================================
!
!BOP
!
! !IROUTINE: gridbox_corners - get coordinates of grid box corners
!
! !INTERFACE:
!
      subroutine gridbox_corners
!
! !DESCRIPTION:
!
! These fields are only used for netcdf history output, and the
! ghost cell values are not needed.
! NOTE:  Extrapolations were used: these fields are approximate!
!
! !REVISION HISTORY:
!
! authors:   A. McLaren, Met Office
!            E. Hunke, LANL
!
! !USES:
      use ice_work, only: work_g2
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
          i,j,iblk,icorner,& ! index counters
          ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      !-------------------------------------------------------------
      ! Get coordinates of grid boxes for each block as follows:
      ! (1) SW corner, (2) SE corner, (3) NE corner, (4) NW corner
      !-------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            latu_bounds(1,i,j,iblk)=TLAT(i  ,j  ,iblk)*rad_to_deg
            latu_bounds(2,i,j,iblk)=TLAT(i+1,j  ,iblk)*rad_to_deg
            latu_bounds(3,i,j,iblk)=TLAT(i+1,j+1,iblk)*rad_to_deg
            latu_bounds(4,i,j,iblk)=TLAT(i  ,j+1,iblk)*rad_to_deg         

            lonu_bounds(1,i,j,iblk)=TLON(i  ,j  ,iblk)*rad_to_deg
            lonu_bounds(2,i,j,iblk)=TLON(i+1,j  ,iblk)*rad_to_deg
            lonu_bounds(3,i,j,iblk)=TLON(i+1,j+1,iblk)*rad_to_deg
            lonu_bounds(4,i,j,iblk)=TLON(i  ,j+1,iblk)*rad_to_deg         

         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !----------------------------------------------------------------
      ! extrapolate on global grid to get edge values
      !----------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global))
      else
         allocate(work_g2(1,1))
      endif

      work1(:,:,:) = latu_bounds(2,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      latu_bounds(2,:,:,:) = work1(:,:,:)

      work1(:,:,:) = latu_bounds(3,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      latu_bounds(3,:,:,:) = work1(:,:,:)

      work1(:,:,:) = latu_bounds(4,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      latu_bounds(4,:,:,:) = work1(:,:,:)

      work1(:,:,:) = lonu_bounds(2,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lonu_bounds(2,:,:,:) = work1(:,:,:)

      work1(:,:,:) = lonu_bounds(3,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lonu_bounds(3,:,:,:) = work1(:,:,:)

      work1(:,:,:) = lonu_bounds(4,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lonu_bounds(4,:,:,:) = work1(:,:,:)

      deallocate(work_g2)

      !----------------------------------------------------------------
      ! Convert longitude to Degrees East >0 for history output
      !----------------------------------------------------------------

      allocate(work_g2(nx_block,ny_block))  ! not used as global here
      !NOTE - this is commented below due to problems with OpenMP
      !reproducibility in this loop  
      !!$OMP PARALLEL DO PRIVATE(iblk,icorner,work_g2)
      do iblk = 1, nblocks
         do icorner = 1, 4
            work_g2(:,:) = lont_bounds(icorner,:,:,iblk) + c360
            where (work_g2 > c360) work_g2 = work_g2 - c360
            where (work_g2 < c0 )  work_g2 = work_g2 + c360
            lont_bounds(icorner,:,:,iblk) = work_g2(:,:)

            work_g2(:,:) = lonu_bounds(icorner,:,:,iblk) + c360
            where (work_g2 > c360) work_g2 = work_g2 - c360
            where (work_g2 < c0 )  work_g2 = work_g2 + c360
            lonu_bounds(icorner,:,:,iblk) = work_g2(:,:)
         enddo
      enddo
      !!$OMP END PARALLEL DO
      deallocate(work_g2)

      end subroutine gridbox_corners

!=======================================================================
!
!BOP
!
! !IROUTINE: gridbox_verts - coordinates of grid box vertices
!
! !INTERFACE:
!
      subroutine gridbox_verts(work_g,vbounds)
!
! !DESCRIPTION:
!
! NOTE:  Boundary conditions for fields on NW, SW, SE corners
!        have not been implemented; using NE corner location for all.
!        Extrapolations are also used: these fields are approximate!
!
! !REVISION HISTORY:
!
! authors:   A. McLaren, Met Office
!            E. Hunke, LANL
!
! !USES:
      use ice_work, only: work_g2
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind), dimension(:,:), intent(in) :: work_g

      real (kind=dbl_kind), &
          dimension(4,nx_block,ny_block,max_blocks), &
          intent(out) :: vbounds

      integer (kind=int_kind) :: &
          i,j,             & ! index counters
          ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global))
      else
         allocate(work_g2(1,1))
      endif

      !-------------------------------------------------------------
      ! Get coordinates of grid boxes for each block as follows:
      ! (1) SW corner, (2) SE corner, (3) NE corner, (4) NW corner
      !-------------------------------------------------------------

      work_g2(:,:) = c0
      if (my_task == master_task) then
         do j = 2, ny_global
         do i = 2, nx_global
            work_g2(i,j) = work_g(i-1,j-1) * rad_to_deg
         enddo
         enddo
         ! extrapolate
         do j = 1, ny_global
            work_g2(1,j) = c2*work_g2(2,j) - work_g2(3,j)
         enddo
         do i = 1, nx_global
            work_g2(i,1) = c2*work_g2(i,2) - work_g2(i,3)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      vbounds(1,:,:,:) = work1(:,:,:)

      work_g2(:,:) = c0
      if (my_task == master_task) then
         do j = 2, ny_global
         do i = 1, nx_global
            work_g2(i,j) = work_g(i,j-1) * rad_to_deg
         enddo
         enddo
         ! extrapolate
         do i = 1, nx_global
            work_g2(i,1) = (c2*work_g2(i,2) - work_g2(i,3))
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      vbounds(2,:,:,:) = work1(:,:,:)

      work_g2(:,:) = c0
      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 1, nx_global
            work_g2(i,j) = work_g(i,j) * rad_to_deg
         enddo
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      vbounds(3,:,:,:) = work1(:,:,:)

      work_g2(:,:) = c0
      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 2, nx_global
            work_g2(i,j) = work_g(i-1,j  ) * rad_to_deg         
         enddo
         enddo
         ! extrapolate
         do j = 1, ny_global
            work_g2(1,j) = c2*work_g2(2,j) - work_g2(3,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      vbounds(4,:,:,:) = work1(:,:,:)

      deallocate (work_g2)

      end subroutine gridbox_verts

!=======================================================================

      end module ice_grid

!=======================================================================
