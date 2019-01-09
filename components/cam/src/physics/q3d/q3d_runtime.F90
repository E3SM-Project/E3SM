module q3d_runtime
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
   private
   save

   ! Public interfaces
   public :: q3d_crm_readnl         ! Read the q3d_nl namelist
   public :: q3d_set_channel_decomp ! Setup channel distribution
   public :: define_crm_grids       ! Set up CRM and vGCM output grids
   public :: crm_chunk              ! Chunk # given a channel # and segment #

   ! Protected values
   character(len=256), public, protected :: q3d_channel_filename    = ''
   integer,            public, protected :: q3d_ncrm                = 0
   logical,            public, protected :: q3d_check_energy_change = .false.
   integer,            public, protected :: q3d_channel_width       = 0
   integer,            public, protected :: q3d_total_channels      = 0
   integer,            public, protected :: q3d_begchan             = 0
   integer,            public, protected :: q3d_endchan             = 0
   integer,            public, protected :: crm_decomp              = 0
   integer,            public, protected :: vGCM_ext_decomp         = 0
   integer,            public, protected :: vGCM_decomp             = 0

   ! Private values
   logical :: ncrm_set   = .false.
   logical :: decomp_set = .false.

CONTAINS

   subroutine q3d_crm_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use spmd_utils,      only: masterproc, masterprocid, mpicom
      use spmd_utils,      only: MPI_CHARACTER, MPI_INTEGER, MPI_LOGICAL
      use cam_logfile,     only: iulog
      use cam_abortutils,  only: endrun

      character(len=*), intent(in) :: nlfile  ! filepath for namelist input file

      ! Local variables
      integer                     :: unitn
      integer                     :: ierr
      integer                     :: ncrm
      integer                     :: q3d_crm_channel_width
      character(len=*), parameter :: subname = 'q3d_crm_readnl'

      namelist /q3d_nl/ q3d_channel_filename, q3d_ncrm,                       &
           q3d_crm_channel_width, q3d_check_energy_change
      !------------------------------------------------------------------------

      ncrm = q3d_ncrm ! Save initial state for cross-check
      q3d_channel_filename = ''
      q3d_ncrm = 0
      q3d_crm_channel_width = 0
      q3d_check_energy_change = .false.
      if (masterproc) then
         write(iulog, *) subname, ": reading q3d_nl namelist..."
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'q3d_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, q3d_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading q3d_nl namelist')
            end if
         else
            ! The q3d_nl is required for now
            call endrun(subname // ':: ERROR, q3d_nl namelist not found')
         end if
         close(unitn)
         call freeunit(unitn)
      end if

      ! Broadcast namelist variables
      call mpi_bcast(q3d_channel_filename, 256, MPI_CHARACTER, masterprocid, mpicom, ierr)
      call mpi_bcast(q3d_ncrm, 1, MPI_INTEGER, masterprocid, mpicom, ierr)
      call mpi_bcast(q3d_check_energy_change, 1, MPI_LOGICAL, masterprocid, mpicom, ierr)
      call mpi_bcast(q3d_crm_channel_width, 1, MPI_INTEGER, masterprocid, mpicom, ierr)

      if (ncrm_set .and. (ncrm /= q3d_ncrm)) then
         call endrun(subname//': q3d_ncrm already set')
      end if
      q3d_channel_width = q3d_crm_channel_width
      ncrm_set = .true.
      if (masterproc) then
         write(iulog, '(2a,i0)') subname, ': q3d_ncrm set to ',q3d_ncrm
         write(iulog, '(2a,i0)') subname, ': q3d_channel_width set to ',q3d_channel_width
         write(iulog, '(3a)') subname, ': Channel info will be read from ',trim(q3d_channel_filename)
         if (q3d_check_energy_change) then
            write(iulog, '(2a)') subname, ': Check that the energy and water change matches boundary fluxes'
         else
            write(iulog, '(2a)') subname, ': No check that the energy and water change matches boundary fluxes'
         end if
      end if

   end subroutine q3d_crm_readnl

   subroutine q3d_set_channel_decomp(total_chans, begchan, endchan)
      use cam_abortutils,  only: endrun

      integer, intent(in) :: total_chans
      integer, intent(in) :: begchan
      integer, intent(in) :: endchan

      if (decomp_set .and. (q3d_total_channels /= total_chans)) then
         call endrun('q3d_set_channel_decomp: incompatible channel decomp')
      end if
      q3d_total_channels = total_chans
      q3d_begchan        = begchan
      q3d_endchan        = endchan

   end subroutine q3d_set_channel_decomp

   subroutine define_crm_grids()
      use physconst,        only: pi
      use cam_abortutils,   only: endrun
      use spmd_utils,       only: masterproc
      use cam_logfile,      only: iulog
      use cam_grid_support, only: cam_grid_register, iMap, new_cam_grid_id
      use cam_grid_support, only: cam_grid_attribute_register
      use cam_grid_support, only: horiz_coord_t, horiz_coord_create
      use parmsld,          only: nVGCM_seg, nVGCM, nhalo_vGCM
      use vGCM_data_types,  only: channel_vGCM

      integer                          :: crm_dim_size
      integer                          :: vGCM_dim_size
      integer                          :: channel_len, chan_seg_len
      integer                          :: local_size
      integer                          :: chan_num, chan_ind, chan_seg
      integer                          :: index, seg_index
      real(r8),            allocatable :: latvals(:)
      real(r8),            allocatable :: lonvals(:)
      integer,             pointer     :: chan_num_vals(:)
      integer,             pointer     :: chan_ind_vals(:)
      integer,             pointer     :: chan_seg_vals(:)
      integer(iMap),       pointer     :: grid_map(:,:)
      integer(iMap),       allocatable :: coord_map(:)
      type(horiz_coord_t), pointer     :: lat_coord
      type(horiz_coord_t), pointer     :: lon_coord
      real(r8),            parameter   :: rad2deg = 180.0_r8 / pi
      character(len=*),    parameter   :: sub = 'DEFINE_CRM_GRIDS'

!#ifdef JUNG_TEST
      integer                          :: halo_len, nn, npo
      real(r8),            pointer     :: hlatvals(:,:)
      real(r8),            pointer     :: hlonvals(:,:)
      character(len=11),   allocatable :: hlatname(:)
      character(len=11),   allocatable :: hlonname(:)
!#endif

      nullify(grid_map)
      nullify(lat_coord)
      nullify(lon_coord)

      ! Are we ready for this?
      if (q3d_total_channels <= 0) then
         call endrun(sub//': ERROR called before channel decomp set')
      end if
      ! Let's be idempotent
      if (crm_decomp /= 0) then
         if (masterproc) then
            write(iulog, *) sub, ': WARNING: CRM grid already defined'
            return
         end if
      end if
      ! For CRM output, we need a column (crm_col) for every CRM cell
      ! For convenient chunking, every channel segment will be a chunk
      ! chunk number given by crm_chunk(chan_num, seg_num)
      channel_len = nVGCM * q3d_ncrm
      chan_seg_len = nVGCM_seg * q3d_ncrm
      crm_dim_size = q3d_total_channels * channel_len
      local_size = (q3d_endchan - q3d_begchan + 1) * channel_len

      allocate(grid_map(3, local_size))
      allocate(coord_map(local_size))
      allocate(lonvals(local_size))
      allocate(latvals(local_size))
      allocate(chan_num_vals(local_size))
      allocate(chan_seg_vals(local_size))
      allocate(chan_ind_vals(local_size))
      index = 0
      do chan_num = q3d_begchan, q3d_endchan
         do chan_ind = 1, channel_len
            chan_seg = ((chan_ind - 1) / chan_seg_len) + 1
            index = index + 1
            coord_map(index) = ((chan_num - 1) * channel_len) + chan_ind
            seg_index = mod((chan_ind - 1), chan_seg_len) + 1
            grid_map(1, index) = seg_index
            grid_map(2, index) = crm_chunk(chan_num, chan_seg)
            grid_map(3, index) = coord_map(index)

            ! We are using the actual lon/lat for CRM cells in a vGCM
            lonvals(index) = channel_vGCM(chan_num)%vGCM_out(chan_seg)%clon(seg_index) * rad2deg
            latvals(index) = channel_vGCM(chan_num)%vGCM_out(chan_seg)%clat(seg_index) * rad2deg
            chan_num_vals(index) = chan_num
            chan_seg_vals(index) = chan_seg
            chan_ind_vals(index) = chan_ind
         end do
      end do
      lon_coord => horiz_coord_create('lon_crm', 'crm_col', crm_dim_size,     &
           'longitude', 'degrees_east', 1, size(lonvals), lonvals,            &
           map=coord_map)
      lat_coord => horiz_coord_create('lat_crm', 'crm_col', crm_dim_size,     &
           'latitude', 'degrees_north', 1, size(latvals), latvals,            &
           map=coord_map)
      deallocate(latvals)
      deallocate(lonvals)
      crm_decomp = new_cam_grid_id()
      call cam_grid_register('crm_grid', crm_decomp, lat_coord, lon_coord,    &
           grid_map, unstruct=.true.)
      ! Now, create attributes for channel number, segment, and index
      call cam_grid_attribute_register('crm_grid', 'channel_number',          &
           'channel number', 'crm_col', chan_num_vals, map=coord_map)
      call cam_grid_attribute_register('crm_grid', 'channel_segment',         &
           'channel segment for each CRM column', 'crm_col', chan_seg_vals,   &
           map=coord_map)
      call cam_grid_attribute_register('crm_grid', 'crm_index',               &
           'channel index for each CRM column', 'crm_col', chan_ind_vals,     &
           map=coord_map)

      ! Cleanup
      deallocate(coord_map)
      nullify(chan_num_vals) ! Belongs to attribute
      nullify(chan_seg_vals) ! Belongs to attribute
      nullify(chan_ind_vals) ! Belongs to attribute
      nullify(grid_map)  ! Belongs to grid
      nullify(lat_coord) ! Belongs to grid
      nullify(lon_coord) ! Belongs to grid

      ! For vGCM output, we need a column (vGCM_col) for every vGCM cell
      ! For convenient chunking, every channel segment will be a chunk
      ! chunk number given by crm_chunk(chan_num, seg_num)
      ! Note that vGCM output includes the halo cells on the ends
      ! Note that while the vGCM arrays start with negative indices, CAM
      !      grid mapping uses lengths so everything below starts at 1

      chan_seg_len = nVGCM_seg + (nhalo_vGCM * 2)
      channel_len = chan_seg_len * 4
      vGCM_dim_size = q3d_total_channels * channel_len
      local_size = (q3d_endchan - q3d_begchan + 1) * channel_len
      allocate(grid_map(3, local_size))
      allocate(coord_map(local_size))
      allocate(lonvals(local_size))
      allocate(latvals(local_size))

      allocate(chan_num_vals(local_size))
      allocate(chan_seg_vals(local_size))
      allocate(chan_ind_vals(local_size))

!#ifdef JUNG_TEST
      halo_len = nhalo_vGCM * 2
      allocate(hlonvals(local_size,halo_len))
      allocate(hlatvals(local_size,halo_len))
      allocate(hlonname(halo_len))
      allocate(hlatname(halo_len))
!#endif

      index = 0
      do chan_num = q3d_begchan, q3d_endchan
         do chan_seg = 1, 4
            do chan_ind = 1, chan_seg_len
               index = index + 1
               coord_map(index) = ((chan_num - 1) * channel_len) + ((chan_seg - 1) * chan_seg_len) + chan_ind
               grid_map(1, index) = chan_ind
               grid_map(2, index) = crm_chunk(chan_num, chan_seg)
               grid_map(3, index) = coord_map(index)
               seg_index = chan_ind - nhalo_vGCM
               lonvals(index) = channel_vGCM(chan_num)%vGCM_state(chan_seg)%lon(0, seg_index) * rad2deg
               latvals(index) = channel_vGCM(chan_num)%vGCM_state(chan_seg)%lat(0, seg_index) * rad2deg

!#ifdef JUNG_TEST
               do nn=1,nhalo_vGCM
                 npo = (nn-1) - nhalo_vGCM
                 hlonvals(index,nn) = &
                         channel_vGCM(chan_num)%vGCM_state(chan_seg)%lon(npo, seg_index) * rad2deg
                 hlatvals(index,nn) = &
                         channel_vGCM(chan_num)%vGCM_state(chan_seg)%lat(npo, seg_index) * rad2deg
               enddo
               do nn=1,nhalo_vGCM
                 npo = nn+nhalo_vGCM
                 hlonvals(index,npo) = &
                         channel_vGCM(chan_num)%vGCM_state(chan_seg)%lon(nn, seg_index) * rad2deg
                 hlatvals(index,npo) = &
                         channel_vGCM(chan_num)%vGCM_state(chan_seg)%lat(nn, seg_index) * rad2deg
               enddo
!#endif
               chan_num_vals(index) = chan_num
               chan_seg_vals(index) = chan_seg
               chan_ind_vals(index) = seg_index
            end do
         end do
      end do
      lon_coord => horiz_coord_create('lon_vGCM_ext', 'vGCM_ext_col',         &
           vGCM_dim_size, 'longitude', 'degrees_east', 1,                     &
           size(lonvals), lonvals, map=coord_map)
      lat_coord => horiz_coord_create('lat_vGCM_ext', 'vGCM_ext_col',         &
           vGCM_dim_size, 'latitude', 'degrees_north', 1,                     &
           size(latvals), latvals, map=coord_map)
      deallocate(latvals)
      deallocate(lonvals)

      vGCM_ext_decomp = new_cam_grid_id()
      call cam_grid_register('vGCM_ext_grid', vGCM_ext_decomp,                &
           lat_coord, lon_coord, grid_map, unstruct=.true.)
      ! Now, create attributes for channel number, segment, and index
      call cam_grid_attribute_register('vGCM_ext_grid',                       &
           'vGCM_e_chan_num', 'vGCM_ext channel number', 'vGCM_ext_col',      &
           chan_num_vals, map=coord_map)
      call cam_grid_attribute_register('vGCM_ext_grid', 'vGCM_e_chan_seg',    &
           'vGCM_ext channel segment for each CRM column', 'vGCM_ext_col',    &
           chan_seg_vals, map=coord_map)
      call cam_grid_attribute_register('vGCM_ext_grid', 'vGCM_e_index',       &
           'vGCM channel index for each CRM column, including ghost cells',   &
           'vGCM_ext_col', chan_ind_vals, map=coord_map)

!#ifdef JUNG_TEST
      do nn=1,halo_len
        write(unit=hlonname(nn),fmt='(a10,i1.1)') 'vGCM_e_lon',nn
        write(unit=hlatname(nn),fmt='(a10,i1.1)') 'vGCM_e_lat',nn

        call cam_grid_attribute_register('vGCM_ext_grid',                &
           hlonname(nn), 'vGCM_ext channel longitude', 'vGCM_ext_col',   &
           hlonvals(:,nn), map=coord_map)

        call cam_grid_attribute_register('vGCM_ext_grid',                &
           hlatname(nn), 'vGCM_ext channel latitude', 'vGCM_ext_col',    &
           hlatvals(:,nn), map=coord_map)
      enddo

      deallocate(hlonname)
      deallocate(hlatname)
      nullify(hlonvals)      ! Belongs to attribute
      nullify(hlatvals)      ! Belongs to attribute
!#endif

      ! Cleanup
      deallocate(coord_map)
      nullify(chan_num_vals) ! Belongs to attribute
      nullify(chan_seg_vals) ! Belongs to attribute
      nullify(chan_ind_vals) ! Belongs to attribute
      nullify(grid_map)  ! Belongs to grid
      nullify(lat_coord) ! Belongs to grid
      nullify(lon_coord) ! Belongs to grid

      ! For vGCM tendency output, we need a column (vGCM_col) for every
      ! vGCM cell.
      ! For convenient chunking, every channel segment will be a chunk
      ! chunk number given by crm_chunk(chan_num, seg_num)
      channel_len = nVGCM
      chan_seg_len = nVGCM_seg
      vGCM_dim_size = q3d_total_channels * channel_len
      local_size = (q3d_endchan - q3d_begchan + 1) * channel_len

      allocate(grid_map(3, local_size))
      allocate(coord_map(local_size))
      allocate(lonvals(local_size))
      allocate(latvals(local_size))
      allocate(chan_num_vals(local_size))
      allocate(chan_seg_vals(local_size))
      allocate(chan_ind_vals(local_size))
      index = 0
      do chan_num = q3d_begchan, q3d_endchan
         do chan_ind = 1, channel_len
            chan_seg = ((chan_ind - 1) / chan_seg_len) + 1
            index = index + 1
            coord_map(index) = ((chan_num - 1) * channel_len) + chan_ind
            seg_index = mod((chan_ind - 1), chan_seg_len) + 1
            grid_map(1, index) = seg_index
            grid_map(2, index) = crm_chunk(chan_num, chan_seg)
            grid_map(3, index) = coord_map(index)
            lonvals(index) = channel_vGCM(chan_num)%vGCM_state(chan_seg)%lon(0, seg_index) * rad2deg
            latvals(index) = channel_vGCM(chan_num)%vGCM_state(chan_seg)%lat(0, seg_index) * rad2deg
            chan_num_vals(index) = chan_num
            chan_seg_vals(index) = chan_seg
            chan_ind_vals(index) = chan_ind
         end do
      end do
      lon_coord => horiz_coord_create('lon_vGCM', 'vGCM_col',                 &
           vGCM_dim_size, 'longitude', 'degrees_east', 1,                     &
           size(lonvals), lonvals, map=coord_map)
      lat_coord => horiz_coord_create('lat_vGCM', 'vGCM_col',                 &
           vGCM_dim_size, 'latitude', 'degrees_north', 1,                     &
           size(latvals), latvals, map=coord_map)
      deallocate(latvals)
      deallocate(lonvals)
      vGCM_decomp = new_cam_grid_id()
      call cam_grid_register('vGCM_grid', vGCM_decomp,                        &
           lat_coord, lon_coord, grid_map, unstruct=.true.)
      ! Now, create attributes for channel number, segment, and index
      call cam_grid_attribute_register('vGCM_grid', 'vGCM_chan_num',          &
           'vGCM channel number', 'vGCM_col', chan_num_vals, map=coord_map)
      call cam_grid_attribute_register('vGCM_grid', 'vGCM_chan_seg',          &
           'vGCM channel segment for each vGCM column', 'vGCM_col',           &
           chan_seg_vals, map=coord_map)
      call cam_grid_attribute_register('vGCM_grid', 'vGCM_index',             &
           'channel index for each vGCM column', 'vGCM_col',                  &
           chan_ind_vals, map=coord_map)

      ! Cleanup
      deallocate(coord_map)
      nullify(chan_num_vals) ! Belongs to attribute
      nullify(chan_seg_vals) ! Belongs to attribute
      nullify(chan_ind_vals) ! Belongs to attribute
      nullify(grid_map)  ! Belongs to grid
      nullify(lat_coord) ! Belongs to grid
      nullify(lon_coord) ! Belongs to grid

   end subroutine define_crm_grids

   integer function crm_chunk(chan_num, seg_num)
      use cam_abortutils,  only: endrun

      integer, intent(in) :: chan_num
      integer, intent(in) :: seg_num

      if ((chan_num < q3d_begchan) .or. (chan_num > q3d_endchan)) then
         call endrun("CRM_CHUNK: chan_num out of bounds")
      end if
      if ((seg_num < 1) .or. (seg_num > 4)) then
         call endrun("CRM_CHUNK: seg_num out of bounds")
      end if
      crm_chunk = ((chan_num - 1) * 4) + seg_num
   end function crm_chunk

end module q3d_runtime
