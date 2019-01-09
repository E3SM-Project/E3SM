module q3d_comm
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use cam_logfile,      only: iulog
  use shr_sys_mod,      only: shr_sys_flush
  use parmsld,          only: nhalo_vgcm, nvgcm_seg, nlevel, ntracer
  use ppgrid,           only: pcols, begchunk, endchunk
  use constituents,     only: pcnst
  use cam_grid_support, only: iMap
  use q3d_coupling,     only: halo_index_t, element_index_t

  implicit none
  private
  save

  !===================
  ! Public interface
  !===================
  public :: q3d_comm_init
  public :: q3d_comm_state_to_vGCM
  public :: q3d_comm_vGCM_to_tend

#ifdef JUNG_TEST
  public :: q3d_comm_vGCM_to_dout
#endif

  ! Derived MPI types
  integer, public, protected :: vGCM_state_mpi_type
  integer, public, protected :: vGCM_tend_mpi_type

  !==============================
  ! Private interfaces and types
  !==============================

  ! Information for each physics (vGCM) column
  type :: vgcm_column_loc_t
     integer(iMap) :: gcid = 0 ! Global column index
     integer  :: task                                                        = 0
     integer  :: face                                                        = 0
     integer  :: segment                                                     = 0
     integer  :: channel1_num(-nhalo_vGCM:nhalo_vGCM)                        = 0
     integer  :: channel1_index                                              = 0
     integer  :: channel2_num(-nhalo_vGCM:nhalo_vGCM)                        = 0
     integer  :: channel2_index                                              = 0
     ! ghost_source contains the global column number assigned to some
     ! of this column's off-face halo points
     integer  :: ghost_source(-nhalo_vGCM:nhalo_vGCM,-nhalo_vGCM:nhalo_vGCM) = 0
     real(r8) :: rlat ! column center latitude (approx for ghosts), radians
     real(r8) :: rlon ! column center longitude (approx for ghosts), radians
  end type vgcm_column_loc_t

  type(vgcm_column_loc_t), allocatable :: phys_cols(:,:) ! pcols, chunks

  ! Information to send from physics to CRM tasks to locate column memory
  type :: phys_to_vgm_loc_t
     integer :: phys_task          ! Where to return vGCM tendencies
     integer :: phys_chunk         ! Used to identify unique column
     integer :: phys_col           ! Used to identify unique column
     integer :: phys_chan_type     ! From column's channel1 or channel 2 info?
     integer :: vgcm_task          ! Where to send vGCM state
     integer :: vgcm_channel_num   ! Channel number on vGCM task
     integer :: vgcm_channel_seg   ! Channel segment number on vGCM task
     integer :: vgcm_channel_index ! Channel index on vGCM task
     integer :: vgcm_ghost_index   ! Ghost index on vGCM task (zero is center)
     integer :: phys_tag           ! msg tag for phys to vgcm
     integer :: vgcm_tag           ! msg tag for vgcm to phys
  end type phys_to_vgm_loc_t

  type(phys_to_vgm_loc_t), allocatable :: phys_sendcols(:)
  integer,                 allocatable :: send_indexing(:)
  integer                              :: num_phys_sendmsgs = 0
  integer                              :: num_vgcm_sendmsgs = 0
  type(phys_to_vgm_loc_t), allocatable :: vGCM_recvcols(:)
  integer,                 allocatable :: recv_indexing(:)
  integer                              :: num_phys_recvmsgs = 0
  integer                              :: num_vgcm_recvmsgs = 0
  integer                              :: nedge = 0 ! # columns per segment
  integer(iMap)                        :: max_prog_column = 0
  ! Which channel segment is each face in (for both channel1 and channel2)?
  integer,                 parameter   :: channel_segment(6, 2) =             &
       reshape((/ 2, 3, 4, 1, 3, 1, 2, 3, 4, 1, 2, 4 /), (/ 6, 2 /))

  type :: phys_to_vgcm_state_t
     real(kind=r8) :: T(nlevel)         ! temperature
     real(kind=r8) :: pint(nlevel+1)    ! dry pressure at interfaces
     real(kind=r8) :: QT(nlevel, pcnst) ! tracer array
     real(kind=r8) :: U(nlevel)         ! zonal velocity component
     real(kind=r8) :: V(nlevel)         ! meridional velocity component
     real(kind=r8) :: omega(nlevel)     ! vertical pressure velocity
     real(kind=r8) :: zm_int(nlevel+1)  ! dry pressure at interfaces
  end type phys_to_vgcm_state_t

  type vgcm_to_phys_tend_t
     real(kind=r8) :: dT(nlevel)         ! temperature tendency
     real(kind=r8) :: dQT(nlevel, pcnst) ! tracer tendencies
     real(kind=r8) :: dU(nlevel)         ! zonal wind tendency
     real(kind=r8) :: dV(nlevel)         ! meridional wind tendency
  end type vgcm_to_phys_tend_t

#ifdef JUNG_TEST
  type vgcm_to_phys_dout_t
     real(kind=r8) :: dT(nlevel)         ! temperature tendency
     real(kind=r8) :: dQT(nlevel, pcnst) ! tracer tendencies
     real(kind=r8) :: dU(nlevel)         ! zonal wind tendency
     real(kind=r8) :: dV(nlevel)         ! meridional wind tendency

     real(kind=r8) :: T(nlevel)          ! temperature
     real(kind=r8) :: QT(nlevel, pcnst)  ! tracers (including water species)
     real(kind=r8) :: U(nlevel)          ! zonal wind
     real(kind=r8) :: V(nlevel)          ! meridional wind
  end type vgcm_to_phys_dout_t
#endif

  ! We need to know the number of columns in all blocks
  integer,       allocatable :: ncols_p(:)        ! Replacement for get_ncols_p
  integer                    :: ghost_chunk_s     ! First chunk of ghost colums
  integer                    :: ghost_chunk_e     ! First chunk of ghost colums
  ! We need to be able to efficiently fill ghost cells during d_p_coupling
  ! ghost_fill_index contains the location of an offset to a off-face halo
  ! column relative to the source (nearest) on-face column
  type(halo_index_t),    allocatable :: ghost_fill_index(:)
  ! ghost_index provides a link from column & chunk to halo transfer index
  integer,               allocatable :: ghost_index(:,:)
  ! We need element information not available from d_p_coupling
  ! element_fill_index contains the dynamics decomposition location for each
  ! on-face column
  type(element_index_t), allocatable :: element_fill_index(:)
  ! element_index provides a link from column & chunk to elem transfer index
  integer,               allocatable :: element_index(:,:)

  ! Fun and games with MPI tag upper bounds
  integer                    :: mpi_tag_max = 32767 ! Min allowed value in MPI
  ! A physically unrealistic value to use in unused slots
  real(r8), parameter        :: badval = HUGE(1.0_r8)

  ! Constituent indices
  integer :: ixcldliq  = -1 ! cloud liquid amount index
  integer :: ixcldice  = -1 ! cloud ice amount index
  integer :: ixrain    = -1 ! rain amount index
  integer :: ixsnow    = -1 ! snow amount index
  integer :: ixgraupel = -1 ! graupel amount index
  ! cnst_indices holds the CAM Q index for each QT member
  integer :: cnst_indices(ntracer) ! Mapping to QT

  interface comm_to_vgcm
     module procedure comm_to_vgcm_int0
     module procedure comm_to_vgcm_r8
     module procedure comm_to_vgcm_state
  end interface comm_to_vgcm

  interface comm_to_phys
     module procedure comm_to_phys_tend
  end interface comm_to_phys

  private :: q3d_comm_read_grid
  private :: q3d_channel_number_to_pe
  private :: get_channel_segment
  private :: cons_phys_tag
  private :: cons_vgcm_tag

  !=======================================================================
contains
  !=======================================================================

  subroutine q3d_comm_init()
    use cam_abortutils,  only: endrun
    use constituents,    only: cnst_get_ind
    use spmd_utils,      only: mpicom, MPI_REAL8, MPI_ADDRESS_KIND, MPI_TAG_UB
    use cam_history_support, only: add_hist_coord, add_vert_coord
    use q3d_runtime,     only: q3d_channel_filename, q3d_begchan, q3d_endchan
    use q3d_coupling,    only: q3d_halo_coord, halo_index_lon, halo_index_lat
    use vGCM_data_types, only: vGCM_allocate_channel_data, channel_vGCM

    ! Local data
    integer(kind=MPI_ADDRESS_KIND)   :: offsets(20)    ! For new MPI types
    integer                          :: origtypes(20)  ! For new MPI types
    integer                          :: lengths(20)    ! For new MPI types
    integer(kind=MPI_ADDRESS_KIND)   :: extent         ! For new MPI types
    integer                          :: num_fields     ! For new MPI types
    logical                          :: flag           ! For MPI query
    integer                          :: ierr           ! Error code
    integer                          :: index, chnk, col, ncol
    integer                          :: chan, seg, vind, gind
    integer                          :: qinds(6)
    integer,  allocatable            :: face_arr(:,:)
    integer,  allocatable            :: vGCM_faces(:)
    real(r8), allocatable            :: fld_coord(:,:)
    real(r8), allocatable            :: coord_arr(:,:)
    real(r8), allocatable            :: vGCM_coords(:)
    integer                          :: dummy_type     ! For new MPI types
    type(phys_to_vgcm_state_t)       :: dummy_state(2) ! For new MPI types
    type(vgcm_to_phys_tend_t)        :: dummy_tend(2)  ! For new MPI types

    ! Find out the correct value for mpi_tag_ub
    call MPI_comm_get_attr(mpicom, MPI_TAG_UB, offsets(1), flag, ierr)
    if (flag) then
       mpi_tag_max = int(offsets(1))
    end if

    ! Retrieve the Q3D constituent indices and figure out QT indices
    index = 1
    qinds(index) = 1 ! Water vapor always 1 in CAM
    call cnst_get_ind('CLDLIQ', ixcldliq)
    index = index + 1
    qinds(index) = ixcldliq
    call cnst_get_ind('CLDICE', ixcldice)
    index = index + 1
    qinds(index) = ixcldice
    call cnst_get_ind('RAINQM', ixrain)
    index = index + 1
    qinds(index) = ixrain
    call cnst_get_ind('SNOWQM', ixsnow)
    index = index + 1
    qinds(index) = ixsnow
    call cnst_get_ind('GRAUPELQM', ixgraupel)
    index = index + 1
    qinds(index) = ixgraupel
    ! cnst_indices holds the CAM Q index for each QT member
    col = 0
    cnst_indices = -1
    do index = 2, pcnst
       if (.not. ANY(qinds == index)) then
          col = col + 1
          if (col > ntracer) then
             call endrun('q3d_comm_init: Not enough room in QT for CAM tracers')
          end if
          cnst_indices(col) = index
       end if
    end do
    if (col < ntracer) then
       call endrun('q3d_comm_init: Too many QT entries for CAM tracers')
    end if

    ! Set up the grid
    call q3d_comm_read_grid(q3d_channel_filename)
    ! Set up the vGCM array
    call vGCM_allocate_channel_data(q3d_begchan, q3d_endchan)

    ! Define the MPI types needed to send state and tendencies between the
    ! vGCM and the standard physics decompositions

    ! type phys_to_vgcm_state_t
    num_fields = 0
    ! All fields are MPI_REAL8
    origtypes(:) = MPI_REAL8
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%T,      offsets(num_fields), ierr)
    lengths(num_fields) = nlevel
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%pint,   offsets(num_fields), ierr)
    lengths(num_fields) = nlevel + 1
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%QT,     offsets(num_fields), ierr)
    lengths(num_fields) = nlevel * pcnst
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%U,      offsets(num_fields), ierr)
    lengths(num_fields) = nlevel
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%V,      offsets(num_fields), ierr)
    lengths(num_fields) = nlevel
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%omega,  offsets(num_fields), ierr)
    lengths(num_fields) = nlevel
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%zm_int, offsets(num_fields), ierr)
    lengths(num_fields) = nlevel + 1
    ! Adjust offsets and create type
    do index = num_fields, 1, -1
       offsets(index) = offsets(index) - offsets(1)
    end do
    call MPI_type_create_struct(num_fields, lengths(1:num_fields),            &
         offsets(1:num_fields), origtypes(1:num_fields), dummy_type, ierr)
    ! Adjust for padding
    call MPI_Get_address(dummy_state(1)%T, offsets(1), ierr)
    call MPI_Get_address(dummy_state(2)%T, offsets(2), ierr)
    extent = offsets(2) - offsets(1)
    call MPI_type_create_resized(dummy_type, 0_MPI_ADDRESS_KIND, extent, &
         vGCM_state_mpi_type, ierr)
    call MPI_type_commit(vGCM_state_mpi_type, ierr)

    ! type vGCM_tend_t
    num_fields = 0
    ! All fields are MPI_REAL8
    origtypes(:) = MPI_REAL8
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dT, offsets(num_fields), ierr)
    lengths(num_fields) = nlevel
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dQT, offsets(num_fields), ierr)
    lengths(num_fields) = nlevel * pcnst
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dU,  offsets(num_fields), ierr)
    lengths(num_fields) = nlevel
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dV,  offsets(num_fields), ierr)
    lengths(num_fields) = nlevel
    do index = num_fields, 1, -1
       offsets(index) = offsets(index) - offsets(1)
    end do
    call MPI_type_create_struct(num_fields, lengths(1:num_fields),            &
         offsets(1:num_fields), origtypes(1:num_fields), dummy_type, ierr)
    ! Adjust for padding
    call MPI_Get_address(dummy_tend(1)%dt, offsets(1), ierr)
    call MPI_Get_address(dummy_tend(2)%dt, offsets(2), ierr)
    extent = offsets(2) - offsets(1)
    call MPI_type_create_resized(dummy_type, 0_MPI_ADDRESS_KIND, extent, &
         vGCM_tend_mpi_type, ierr)
    call MPI_type_commit(vGCM_tend_mpi_type, ierr)

    !! Send static info to vGCM
    ! Face info
    allocate(vGCM_faces(size(vGCM_recvcols)))
    allocate(face_arr(pcols, begchunk:ghost_chunk_e))
    do chnk = begchunk, ghost_chunk_e
       ncol = ncols_p(chnk)
       do col = 1, ncol
          face_arr(col, chnk) = phys_cols(col, chnk)%face
       end do
    end do
    call comm_to_vgcm(face_arr, vGCM_faces)
    col = 0
    do index = 1, size(vGCM_recvcols)
       ! First recvcol should have phys_tag >= 0 (checked earlier)
       if (vGCM_recvcols(index)%phys_tag >= 0) then
          col = col + 1
       end if
       chan = vGCM_recvcols(index)%vgcm_channel_num
       seg = vGCM_recvcols(index)%vgcm_channel_seg
       channel_vGCM(chan)%vGCM_state(seg)%nface = vGCM_faces(col)
    end do
    ! Cleanup
    deallocate(vGCM_faces)
    deallocate(face_arr)
    ! lat/lon
    allocate(fld_coord(2, size(element_fill_index) + size(ghost_fill_index)))
    allocate(vGCM_coords(size(vGCM_recvcols)))
    allocate(coord_arr(pcols, begchunk:ghost_chunk_e))
    call q3d_halo_coord(element_fill_index, ghost_fill_index, fld_coord)
    do chnk = begchunk, ghost_chunk_e
       ncol = ncols_p(chnk)
       do col = 1, ncol
          if (chnk > endchunk) then
             index = ghost_index(col, chnk) + size(element_fill_index)
          else
             index = element_index(col, chnk)
          end if
          coord_arr(col, chnk) = fld_coord(halo_index_lat, index)
       end do
    end do
    vGCM_coords(:) = badval
    call comm_to_vgcm(coord_arr, vGCM_coords)
    col = 0
    do index = 1, size(vGCM_recvcols)
       ! First recvcol should have phys_tag >= 0 (checked earlier)
       if (vGCM_recvcols(index)%phys_tag >= 0) then
          col = col + 1
       end if
       chan = vGCM_recvcols(index)%vgcm_channel_num
       seg = vGCM_recvcols(index)%vgcm_channel_seg
       vind = vGCM_recvcols(index)%vgcm_channel_index
       gind = vGCM_recvcols(index)%vgcm_ghost_index
       channel_vGCM(chan)%vGCM_state(seg)%lat(gind,vind) = vGCM_coords(col)
    end do
    do chnk = begchunk, ghost_chunk_e
       ncol = ncols_p(chnk)
       do col = 1, ncol
          if (chnk > endchunk) then
             index = ghost_index(col, chnk) + size(element_fill_index)
          else
             index = element_index(col, chnk)
          end if
          coord_arr(col, chnk) = fld_coord(halo_index_lon, index)
       end do
    end do
    vGCM_coords(:) = badval
    call comm_to_vgcm(coord_arr, vGCM_coords)
    col = 0
    do index = 1, size(vGCM_recvcols)
       ! First recvcol should have phys_tag >= 0 (checked earlier)
       if (vGCM_recvcols(index)%phys_tag >= 0) then
          col = col + 1
       end if
       chan = vGCM_recvcols(index)%vgcm_channel_num
       seg = vGCM_recvcols(index)%vgcm_channel_seg
       vind = vGCM_recvcols(index)%vgcm_channel_index
       gind = vGCM_recvcols(index)%vgcm_ghost_index
       channel_vGCM(chan)%vGCM_state(seg)%lon(gind,vind) = vGCM_coords(col)
    end do
    ! Cleanup
    deallocate(vGCM_coords)
    deallocate(coord_arr)

  end subroutine q3d_comm_init

  subroutine get_z_interface(ncol, nlev, pint, Temp, qv, phis, zint)
    use physconst,      only: rair, epsilo, gravit
    use shr_vmath_mod,  only: shr_vmath_log
    use cam_abortutils, only: endrun

    integer,  intent(in)  :: ncol               ! Horizontal array extent
    integer,  intent(in)  :: nlev               ! number of levels
    real(r8), intent(in)  :: pint(nlev+1, ncol) ! full pressure at interfaces
    real(r8), intent(in)  :: qv  (nlev,   ncol) ! dry water vapor mixing ratio
    real(r8), intent(in)  :: Temp(nlev,   ncol) ! temperature
    real(r8), intent(in)  :: PHIS(ncol)         ! surface geopotential

    real(r8), intent(out) :: zint(nlev+1, ncol) ! height at interface
    !
    !---------------------------Local variables-----------------------------
    !
    integer  :: i, k
    real(r8) :: hkl(ncol)              ! off-diagonal element
    real(r8) :: rog                    ! Rair / gravit
    real(r8) :: tv(nlev, ncol)         ! Virtual temperature
    real(r8) :: log_pint(nlev+1, ncol) ! log of pint
    real(r8) :: inv_epsilon

    rog = rair / gravit
    inv_epsilon = 1.0_r8 / Epsilo
    !
    ! compute virtual temperature
    !
    do i = 1, ncol
       do k = 1, nlev
          tv(k,i) = Temp(k,i) * (1.0_r8 + inv_epsilon * qv(k,i)) / (1.0_r8 + qv(k,i))
       end do
    end do
    !
    ! compute log of pint
    !
    if (ANY(pint < 1.0e-12_r8)) then
       call endrun('get_z_interface: Bad pint value(s)')
    end if
    call shr_vmath_log(pint, log_pint, ncol*(nlev+1))
    !
    ! Compute zi from bottom up.
    !
    zint(nlev+1, :) = phis(:) / gravit
    do k = nlev, 1, -1
      hkl(:) = (log_pint(k+1, :) - log_pint(k, :))
      ! dynamics uses hkl(:,:) = pdel(:,:,k) / pmid(:,:,k)
      zint(k, :) = zint(k+1, :) + (rog * tv(k, :) * hkl(:))
    end do
  end subroutine get_z_interface

  subroutine q3d_comm_state_to_vGCM(phys_state)
    use cam_abortutils,  only: endrun
    use ppgrid,          only: pver, pverp
    use physics_types,   only: physics_state
    use q3d_coupling,    only: q3d_halo_data, halo_num_surface
    use q3d_coupling,    only: halo_index_u_con, halo_index_v_con
    use q3d_coupling,    only: halo_index_Temp, halo_index_OMEGA
    use q3d_coupling,    only: halo_index_Q, halo_num_midpoint
    use q3d_coupling,    only: halo_index_pint, halo_index_zi
    use q3d_coupling,    only: halo_num_interface, halo_index_PHIS
    use q3d_runtime,     only: q3d_check_energy_change, q3d_begchan, q3d_endchan
    use vGCM_data_types, only: channel_vGCM

    type(physics_state), intent(in)         :: phys_state(begchunk:endchunk)

    integer                                 :: chnk, col, ncol, index, mind
    integer                                 :: chan, seg, ind, ghind
    real(kind=r8),              allocatable :: fld_phys_s(:,:)
    real(kind=r8),              allocatable :: fld_phys_m(:,:,:)
    real(kind=r8),              allocatable :: fld_phys_i(:,:,:)
    type(phys_to_vgcm_state_t), allocatable :: phys_sarr(:,:)
    type(phys_to_vgcm_state_t), allocatable :: vGCM_state(:)
    character(len=128)                      :: errmsg
    character(len=*),           parameter   :: sub = 'q3d_comm_state_to_vGCM'

    ! Gather up the halo points for the ghost columns
    ! The first block belongs to prognostic columns followed by ghost columns
    index = size(element_fill_index) + size(ghost_fill_index)
    ! the vGCM may have fewer levels than the dycore but not more
    if (nlevel > pver) then
       call endrun(sub//': nlevel must be <= pver')
    end if
    allocate(fld_phys_s(halo_num_surface, index))
    allocate(fld_phys_m(pver, halo_num_midpoint, index))
    allocate(fld_phys_i(pverp, halo_num_interface, index))

    allocate(vGCM_state(size(vGCM_recvcols)))
    allocate(phys_sarr(pcols, begchunk:ghost_chunk_e))
    ! Grab vGCM data from dynamics including ghost columns
    call q3d_halo_data(element_fill_index, ghost_fill_index, fld_phys_s, fld_phys_m, fld_phys_i)

    ! Possibly check to add the heating rate required for global mean
    if (q3d_check_energy_change) then
       call endrun(sub//': Q3D check energy change not implemented')
    end if

    ! Calculate interface geopotential height
    index = size(element_fill_index) + size(ghost_fill_index)
    call get_z_interface(index, nlevel,                                       &
         fld_phys_i(1:nlevel+1, halo_index_pint, :),                          &
         fld_phys_m(1:nlevel,   halo_index_Temp, :),                          &
         fld_phys_m(1:nlevel,   halo_index_Q,    :),                          &
         fld_phys_s(            halo_index_PHIS, :),                          &
         fld_phys_i(1:nlevel+1, halo_index_zi,   :))

    ! Fill in prognostic data
    do chnk = begchunk, endchunk
       do col = 1, ncols_p(chnk)
          index = element_index(col, chnk)
          phys_sarr(col,chnk)%T(1:nlevel)      = fld_phys_m(1:nlevel,halo_index_Temp,index)
          phys_sarr(col,chnk)%U(1:nlevel)      = fld_phys_m(1:nlevel,halo_index_u_con,index)
          phys_sarr(col,chnk)%V(1:nlevel)      = fld_phys_m(1:nlevel,halo_index_v_con,index)
          phys_sarr(col,chnk)%omega(1:nlevel)  = fld_phys_m(1:nlevel,halo_index_OMEGA,index)
          do mind = 1, pcnst
             phys_sarr(col,chnk)%QT(1:nlevel, mind) = fld_phys_m(1:nlevel,halo_index_Q+mind-1,index)
          end do
          phys_sarr(col,chnk)%pint(1:nlevel+1) = fld_phys_i(1:nlevel+1,halo_index_pint,index)
          phys_sarr(col,chnk)%zm_int(1:nlevel+1) = fld_phys_i(1:nlevel+1,halo_index_zi,index)
       end do
    end do
    ! Now, fill ghost column data
    do chnk = ghost_chunk_s, ghost_chunk_e
       do col = 1, ncols_p(chnk)
          ind = ghost_index(col, chnk)
          index = ind + size(element_fill_index)
          phys_sarr(col,chnk)%T(1:nlevel)      = fld_phys_m(1:nlevel,halo_index_Temp,index)
          phys_sarr(col,chnk)%U(1:nlevel)      = fld_phys_m(1:nlevel,halo_index_u_con,index)
          phys_sarr(col,chnk)%V(1:nlevel)      = fld_phys_m(1:nlevel,halo_index_v_con,index)
          phys_sarr(col,chnk)%omega(1:nlevel)  = fld_phys_m(1:nlevel,halo_index_OMEGA,index)
          do mind = 1, pcnst
             phys_sarr(col,chnk)%QT(1:nlevel, mind) = fld_phys_m(1:nlevel,halo_index_Q+mind-1,index)
          end do
          phys_sarr(col,chnk)%pint(1:nlevel+1) = fld_phys_i(1:nlevel+1,halo_index_pint,index)
          phys_sarr(col,chnk)%zm_int(1:nlevel+1) = fld_phys_i(1:nlevel+1,halo_index_zi,index)
       end do
    end do
    ! Send data from the dynamics / vGCM decomp to the CRM decomp
    call comm_to_vgcm(phys_sarr, vGCM_state)
    ! Initialize unused data slots to zero
    do chan = q3d_begchan, q3d_endchan
       do seg = 1, 4
          do ghind = -nhalo_vGCM, nhalo_vGCM, 2*nhalo_vGCM
             do ind = 1-nhalo_vGCM, nVGCM_seg+nhalo_vGCM, nVGCM_seg+(2*nhalo_vGCM)-1
                channel_vGCM(chan)%vGCM_state(seg)%T(ghind,ind,:) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%pint(ghind,ind,:+1) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%zm_int(ghind,ind,:+1) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%U(ghind,ind,:) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%V(ghind,ind,:) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%omega(ghind,ind,:) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%QV(ghind,ind,:) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%QC(ghind,ind,:) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%QI(ghind,ind,:) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%QR(ghind,ind,:) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%QS(ghind,ind,:) = 0.0_r8
                channel_vGCM(chan)%vGCM_state(seg)%QG(ghind,ind,:) = 0.0_r8
                do mind = 1, ntracer
                   channel_vGCM(chan)%vGCM_state(seg)%QT(ghind,ind,:,mind) = 0.0_r8
                end do
             end do
          end do
       end do
    end do
    ! Unpack the real data into the vGCM structures
    col = 0
    do index = 1, size(vGCM_recvcols)
       ! First recvcol should have phys_tag >= 0 (checked earlier)
       if (vGCM_recvcols(index)%phys_tag >= 0) then
          col = col + 1
       end if
       chan = vGCM_recvcols(index)%vgcm_channel_num
       seg = vGCM_recvcols(index)%vgcm_channel_seg
       ind = vGCM_recvcols(index)%vgcm_channel_index
       ghind = vGCM_recvcols(index)%vgcm_ghost_index
       channel_vGCM(chan)%vGCM_state(seg)%T(ghind,ind,1:nlevel) = vGCM_state(col)%T(1:nlevel)
       channel_vGCM(chan)%vGCM_state(seg)%pint(ghind,ind,1:nlevel+1) = vGCM_state(col)%pint(1:nlevel+1)
       channel_vGCM(chan)%vGCM_state(seg)%zm_int(ghind,ind,1:nlevel+1) = vGCM_state(col)%zm_int(1:nlevel+1)
       channel_vGCM(chan)%vGCM_state(seg)%U(ghind,ind,1:nlevel) = vGCM_state(col)%U(1:nlevel)
       channel_vGCM(chan)%vGCM_state(seg)%V(ghind,ind,1:nlevel) = vGCM_state(col)%V(1:nlevel)
       channel_vGCM(chan)%vGCM_state(seg)%omega(ghind,ind,1:nlevel) = vGCM_state(col)%omega(1:nlevel)
       channel_vGCM(chan)%vGCM_state(seg)%QV(ghind,ind,1:nlevel) = vGCM_state(col)%QT(1:nlevel, 1)
       channel_vGCM(chan)%vGCM_state(seg)%QC(ghind,ind,1:nlevel) = vGCM_state(col)%QT(1:nlevel,ixcldliq)
       channel_vGCM(chan)%vGCM_state(seg)%QI(ghind,ind,1:nlevel) = vGCM_state(col)%QT(1:nlevel,ixcldice)
       channel_vGCM(chan)%vGCM_state(seg)%QR(ghind,ind,1:nlevel) = vGCM_state(col)%QT(1:nlevel,ixrain)
       channel_vGCM(chan)%vGCM_state(seg)%QS(ghind,ind,1:nlevel) = vGCM_state(col)%QT(1:nlevel,ixsnow)
       channel_vGCM(chan)%vGCM_state(seg)%QG(ghind,ind,1:nlevel) = vGCM_state(col)%QT(1:nlevel,ixgraupel)
       do mind = 1, ntracer
          channel_vGCM(chan)%vGCM_state(seg)%QT(ghind,ind,1:nlevel,mind) = vGCM_state(col)%QT(1:nlevel,cnst_indices(mind))
       end do
    end do
    ! Cleanup
    deallocate(fld_phys_s)
    deallocate(fld_phys_m)
    deallocate(fld_phys_i)
    deallocate(vGCM_state)
    deallocate(phys_sarr)

  end subroutine q3d_comm_state_to_vGCM

!#ifdef JUNG_TEST
!  subroutine q3d_comm_vGCM_to_tend(phys_state, ztodt, phys_tend, phys_dout)
!#endif
  subroutine q3d_comm_vGCM_to_tend(phys_state, ztodt, phys_tend)

!#ifdef JUNG_TEST
!    use physics_types,   only: physics_state, physics_tend, physics_dout
!#endif
    use physics_types,   only: physics_state, physics_tend
    use vGCM_data_types, only: channel_vGCM
    use cam_abortutils,  only: endrun

    real(r8), intent(in) :: ztodt  ! physics time step unless nstep=0
    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(physics_tend ), intent(inout) :: phys_tend(begchunk:endchunk)

#ifdef JUNG_TEST
    type(physics_dout),  intent(inout) :: phys_dout(begchunk:endchunk)
#endif

    integer                                 :: chnk, col, ncol, index, sindex
    integer                                 :: chan, seg, cind, mind
    real(r8)                                :: qval, temp(nlevel)
    type(vgcm_to_phys_tend_t),  allocatable :: phys_tarr_c1(:,:)
    type(vgcm_to_phys_tend_t),  allocatable :: phys_tarr_c2(:,:)
    type(vgcm_to_phys_tend_t),  allocatable :: vGCM_tend(:)
    character(len=128)                      :: errmsg
    character(len=*),           parameter   :: sub = 'q3d_comm_vGCM_to_tend'
    
    LOGICAL :: WATER_TEND = .TRUE.      ! JUNG_DEBUG

    allocate(vGCM_tend(num_vgcm_sendmsgs))
    allocate(phys_tarr_c1(pcols, begchunk:endchunk))
    allocate(phys_tarr_c2(pcols, begchunk:endchunk))
    ! Load the tendency data into vGCM_tend
    sindex = 0
    do index = 1, size(vGCM_recvcols)
       if (vGCM_recvcols(index)%vgcm_tag >= 0) then
          sindex = sindex + 1
          chan = vGCM_recvcols(index)%vgcm_channel_num
          seg = vGCM_recvcols(index)%vgcm_channel_seg
          cind = vGCM_recvcols(index)%vgcm_channel_index
          vGCM_tend(sindex)%dT(1:nlevel)           = channel_vGCM(chan)%vGCM_tend(seg)%dT(cind,1:nlevel)
          vGCM_tend(sindex)%dU(1:nlevel)           = channel_vGCM(chan)%vGCM_tend(seg)%dU(cind,1:nlevel)
          vGCM_tend(sindex)%dV(1:nlevel)           = channel_vGCM(chan)%vGCM_tend(seg)%dV(cind,1:nlevel)
          vGCM_tend(sindex)%dQT(1:nlevel,1)         = channel_vGCM(chan)%vGCM_tend(seg)%dQV(cind,1:nlevel)
          vGCM_tend(sindex)%dQT(1:nlevel,ixcldliq)  = channel_vGCM(chan)%vGCM_tend(seg)%dQC(cind,1:nlevel)
          vGCM_tend(sindex)%dQT(1:nlevel,ixcldice)  = channel_vGCM(chan)%vGCM_tend(seg)%dQI(cind,1:nlevel)
          vGCM_tend(sindex)%dQT(1:nlevel,ixrain)    = channel_vGCM(chan)%vGCM_tend(seg)%dQR(cind,1:nlevel)
          vGCM_tend(sindex)%dQT(1:nlevel,ixsnow)    = channel_vGCM(chan)%vGCM_tend(seg)%dQS(cind,1:nlevel)
          vGCM_tend(sindex)%dQT(1:nlevel,ixgraupel) = channel_vGCM(chan)%vGCM_tend(seg)%dQG(cind,1:nlevel)
          do mind = 1, ntracer
             vGCM_tend(sindex)%dQT(1:nlevel,cnst_indices(mind)) = channel_vGCM(chan)%vGCM_tend(seg)%dQT(cind,1:nlevel,mind)
          end do
       end if
    end do
    call comm_to_phys(vGCM_tend, phys_tarr_c1, phys_tarr_c2)
    do chnk = begchunk, endchunk
       ncol = phys_state(chnk)%ncol
       if (ncol /= ncols_p(chnk)) then
          write(errmsg, '(3(a,i0))') ': ncol error for ',chnk,', ',ncol,' /= ',ncols_p(chnk)
          call endrun(sub//trim(errmsg))
       end if
       do col = 1, ncol
          temp(1:nlevel) = (phys_tarr_c1(col, chnk)%dT(1:nlevel) + phys_tarr_c2(col, chnk)%dT(1:nlevel)) / 2.0_r8
          phys_tend(chnk)%dTdt(col,1:nlevel) = temp(1:nlevel)

#ifdef JUNG_TEST
          phys_dout(chnk)%dT(col,1:nlevel)   = temp(1:nlevel)
#endif

          temp(1:nlevel) = (phys_tarr_c1(col, chnk)%dU(1:nlevel) + phys_tarr_c2(col, chnk)%dU(1:nlevel)) / 2.0_r8
          phys_tend(chnk)%dUdt(col,1:nlevel) = temp(1:nlevel)

#ifdef JUNG_TEST
          phys_dout(chnk)%dU(col,1:nlevel)   = temp(1:nlevel)
#endif

          temp(1:nlevel) = (phys_tarr_c1(col, chnk)%dV(1:nlevel) + phys_tarr_c2(col, chnk)%dV(1:nlevel)) / 2.0_r8
          phys_tend(chnk)%dVdt(col,1:nlevel) = temp(1:nlevel)

#ifdef JUNG_TEST
          phys_dout(chnk)%dV(col,1:nlevel)   = temp(1:nlevel)
#endif

      IF (WATER_TEND) THEN
      
          do mind = 1, pcnst
             do index = 1, nlevel
                ! Add Q tendency to state
                qval = phys_state(chnk)%Q(col, index, mind)
                temp(1) = (phys_tarr_c1(col, chnk)%dQT(index,mind) + phys_tarr_c2(col, chnk)%dQT(index,mind)) / 2.0_r8
                phys_state(chnk)%Q(col,index,mind) = qval + temp(1)*ztodt
                
#ifdef JUNG_TEST
                phys_dout(chnk)%dQT(col,index,mind)= temp(1)
#endif
             end do
          end do
      
      ELSE ! WATER_TEND  
              
          ! qv: no change         
          do mind = 1, 1
             do index = 1, nlevel
                ! Add Q tendency to state
                qval = phys_state(chnk)%Q(col, index, mind)
                temp(1) = (phys_tarr_c1(col, chnk)%dQT(index,mind) + phys_tarr_c2(col, chnk)%dQT(index,mind)) / 2.0_r8
                
                phys_state(chnk)%Q(col,index,mind) = qval + temp(1)*ztodt

#ifdef JUNG_TEST
                phys_dout(chnk)%dQT(col,index,mind)= temp(1)
#endif
             end do
          end do   
          
          ! qc, qi, qr, qs, qg: use the updated values
          do mind = 2, pcnst
             do index = 1, nlevel
                ! Add Q tendency to state
                qval = phys_state(chnk)%Q(col, index, mind)
                temp(1) = (phys_tarr_c1(col, chnk)%dQT(index,mind) + phys_tarr_c2(col, chnk)%dQT(index,mind)) / 2.0_r8
                
                phys_state(chnk)%Q(col,index,mind) = temp(1)

#ifdef JUNG_TEST
                phys_dout(chnk)%dQT(col,index,mind)= temp(1)
#endif
             end do
          end do             
          
      ENDIF   ! WATER_TEND            

! JUNG_TEST (positive definite)
          ! water vapor
          do mind = 1, 1
             do index = 1, nlevel
                phys_state(chnk)%Q(col,index,mind) = max(1.e-12_r8,phys_state(chnk)%Q(col,index,mind))
             end do
          end do
          ! others
          do mind = 2, pcnst
             do index = 1, nlevel
                phys_state(chnk)%Q(col,index,mind) = max(0.0_r8,phys_state(chnk)%Q(col,index,mind))
             end do
          end do
! JUNG_TEST          
          
       end do
    end do
  end subroutine q3d_comm_vGCM_to_tend

#ifdef JUNG_TEST

  subroutine q3d_comm_vGCM_to_dout(phys_state, phys_dout)

    use physics_types,   only: physics_state,physics_dout
    use vGCM_data_types, only: channel_vGCM
    use cam_abortutils,  only: endrun

    type(physics_state),   intent(in) :: phys_state(begchunk:endchunk)
    type(physics_dout), intent(inout) :: phys_dout(begchunk:endchunk)

    integer                                 :: chnk, col, ncol, index, sindex
    integer                                 :: chan, seg, cind, mind
    real(r8)                                :: temp(nlevel)
    type(vgcm_to_phys_dout_t),  allocatable :: phys_tarr_c1(:,:)
    type(vgcm_to_phys_dout_t),  allocatable :: phys_tarr_c2(:,:)
    type(vgcm_to_phys_dout_t),  allocatable :: vGCM_out(:)
    character(len=128)                      :: errmsg
    character(len=*),           parameter   :: sub = 'q3d_comm_vGCM_to_dout'

    allocate(vGCM_out(num_vgcm_sendmsgs))
    allocate(phys_tarr_c1(pcols, begchunk:endchunk))
    allocate(phys_tarr_c2(pcols, begchunk:endchunk))
    ! Load the diagnostic data into vGCM_out
    sindex = 0
    do index = 1, size(vGCM_recvcols)
       if (vGCM_recvcols(index)%vgcm_tag >= 0) then
          sindex = sindex + 1
          chan = vGCM_recvcols(index)%vgcm_channel_num
          seg = vGCM_recvcols(index)%vgcm_channel_seg
          cind = vGCM_recvcols(index)%vgcm_channel_index
          vGCM_out(sindex)%T(1:nlevel)            = channel_vGCM(chan)%vGCM_out(seg)%T(cind,1:nlevel)
          vGCM_out(sindex)%U(1:nlevel)            = channel_vGCM(chan)%vGCM_out(seg)%U(cind,1:nlevel)
          vGCM_out(sindex)%V(1:nlevel)            = channel_vGCM(chan)%vGCM_out(seg)%V(cind,1:nlevel)
          vGCM_out(sindex)%QT(1:nlevel,1)         = channel_vGCM(chan)%vGCM_out(seg)%QV(cind,1:nlevel)
          vGCM_out(sindex)%QT(1:nlevel,ixcldliq)  = channel_vGCM(chan)%vGCM_out(seg)%QC(cind,1:nlevel)
          vGCM_out(sindex)%QT(1:nlevel,ixcldice)  = channel_vGCM(chan)%vGCM_out(seg)%QI(cind,1:nlevel)
          vGCM_out(sindex)%QT(1:nlevel,ixrain)    = channel_vGCM(chan)%vGCM_out(seg)%QR(cind,1:nlevel)
          vGCM_out(sindex)%QT(1:nlevel,ixsnow)    = channel_vGCM(chan)%vGCM_out(seg)%QS(cind,1:nlevel)
          vGCM_out(sindex)%QT(1:nlevel,ixgraupel) = channel_vGCM(chan)%vGCM_out(seg)%QG(cind,1:nlevel)
          do mind = 1, ntracer
             vGCM_out(sindex)%QT(1:nlevel,cnst_indices(mind)) = channel_vGCM(chan)%vGCM_out(seg)%QT(cind,1:nlevel,mind)
          end do
       end if
    end do
    call comm_to_phys(vGCM_out, phys_tarr_c1, phys_tarr_c2)
    do chnk = begchunk, endchunk
       ncol = phys_state(chnk)%ncol
       if (ncol /= ncols_p(chnk)) then
          write(errmsg, '(3(a,i0))') ': ncol error for ',chnk,', ',ncol,' /= ',ncols_p(chnk)
          call endrun(sub//trim(errmsg))
       end if
       do col = 1, ncol
          temp(1:nlevel) = (phys_tarr_c1(col, chnk)%T(1:nlevel) + phys_tarr_c2(col, chnk)%T(1:nlevel)) / 2.0_r8
          phys_dout(chnk)%T(col,1:nlevel) = temp(1:nlevel)

          temp(1:nlevel) = (phys_tarr_c1(col, chnk)%U(1:nlevel) + phys_tarr_c2(col, chnk)%U(1:nlevel)) / 2.0_r8
          phys_dout(chnk)%U(col,1:nlevel) = temp(1:nlevel)

          temp(1:nlevel) = (phys_tarr_c1(col, chnk)%V(1:nlevel) + phys_tarr_c2(col, chnk)%V(1:nlevel)) / 2.0_r8
          phys_dout(chnk)%V(col,1:nlevel) = temp(1:nlevel)
          do mind = 1, pcnst
            temp(1:nlevel) = (phys_tarr_c1(col, chnk)%QT(1:nlevel,mind) + phys_tarr_c2(col, chnk)%QT(1:nlevel,mind)) / 2.0_r8
            phys_dout(chnk)%QT(col,1:nlevel,mind) = temp(1:nlevel)
          end do
       end do
    end do
  end subroutine q3d_comm_vGCM_to_dout

#endif

  !====================================================
  ! Private routines
  !====================================================

  subroutine comm_to_vgcm_int0(data_in, data_out)
    use cam_abortutils,  only: endrun
    use spmd_utils,      only: mpicom, MPI_INTEGER, MPI_STATUS_SIZE

    integer, intent(in)  :: data_in(:,:)
    integer, intent(out) :: data_out(:)

    integer :: lbnd_in1, ubnd_in1, lbnd_in2, ubnd_in2, lbnd_out, ubnd_out
    integer :: col, chnk, index
    integer :: choff, cooff ! Offsets in case lower bounds do not match
    integer :: msg_num
    integer :: ierr
    integer :: send_request(num_phys_sendmsgs)
    integer :: recv_request(num_vgcm_recvmsgs)
    integer :: recv_status(MPI_STATUS_SIZE, num_vgcm_recvmsgs)
    character(len=*), parameter :: sub = "comm_to_vgcm_int0"

    ! Make sure the dummy arrays match expectations
    lbnd_in1 = LBOUND(data_in, 1)
    ubnd_in1 = UBOUND(data_in, 1)
    lbnd_in2 = LBOUND(data_in, 2)
    ubnd_in2 = UBOUND(data_in, 2)
    lbnd_out = LBOUND(data_out, 1)
    ubnd_out = UBOUND(data_out, 1)
    if ((ubnd_in1 - lbnd_in1 + 1) /= pcols) then
       call endrun(sub//": first index of data_in wrong")
    end if
    cooff = lbnd_in1 - 1
    if ((ubnd_in2 - lbnd_in2) /= (ghost_chunk_e - begchunk)) then
       call endrun(sub//": second index of data_in wrong")
    end if
    choff = lbnd_in2 - begchunk
    if ((ubnd_out - lbnd_out + 1) /= size(vGCM_recvcols)) then
       call endrun(sub//": size of data_out wrong")
    end if

    ! Post receives
    msg_num = 0
    do index = 1, size(vGCM_recvcols)
       if (vGCM_recvcols(index)%phys_tag >= 0) then
          msg_num = msg_num + 1
          call MPI_irecv(data_out(lbnd_out + msg_num - 1), 1, MPI_INTEGER,    &
               vGCM_recvcols(index)%phys_task, vGCM_recvcols(index)%phys_tag, &
               mpicom, recv_request(msg_num), ierr)
       end if
    end do
    ! Post sends
    msg_num = 0
    do index = 1, size(phys_sendcols)
       if (phys_sendcols(index)%phys_tag >= 0) then
          chnk = phys_sendcols(index)%phys_chunk
          col = phys_sendcols(index)%phys_col
          msg_num = msg_num + 1
          call MPI_isend(data_in(col+cooff, chnk+choff), 1, MPI_INTEGER,      &
               phys_sendcols(index)%vgcm_task,                                &
               phys_sendcols(index)%phys_tag,                                 &
               mpicom, send_request(msg_num), ierr)
       end if
    end do

    ! Process received messages
    do index = 1, num_vgcm_recvmsgs
       call MPI_waitall(num_vgcm_recvmsgs, recv_request, recv_status, ierr)
    end do

  end subroutine comm_to_vgcm_int0

  subroutine comm_to_vgcm_r8(data_in, data_out)
    use cam_abortutils,  only: endrun
    use spmd_utils,      only: mpicom, MPI_REAL8, MPI_STATUS_SIZE

    real(r8), intent(in)  :: data_in(:,:)
    real(r8), intent(out) :: data_out(:)

    integer :: lbnd_in1, ubnd_in1, lbnd_in2, ubnd_in2, lbnd_out, ubnd_out
    integer :: col, chnk, index
    integer :: choff, cooff ! Offsets in case lower bounds do not match
    integer :: msg_num
    integer :: ierr
    integer :: send_request(num_phys_sendmsgs)
    integer :: recv_request(num_vgcm_recvmsgs)
    integer :: recv_status(MPI_STATUS_SIZE, num_vgcm_recvmsgs)
    character(len=*), parameter :: sub = "comm_to_vgcm_r8"

    ! Make sure the dummy arrays match expectations
    lbnd_in1 = LBOUND(data_in, 1)
    ubnd_in1 = UBOUND(data_in, 1)
    lbnd_in2 = LBOUND(data_in, 2)
    ubnd_in2 = UBOUND(data_in, 2)
    lbnd_out = LBOUND(data_out, 1)
    ubnd_out = UBOUND(data_out, 1)
    if ((ubnd_in1 - lbnd_in1 + 1) /= pcols) then
       call endrun(sub//": first index of data_in wrong")
    end if
    cooff = lbnd_in1 - 1
    if ((ubnd_in2 - lbnd_in2) /= (ghost_chunk_e - begchunk)) then
       call endrun(sub//": second index of data_in wrong")
    end if
    choff = lbnd_in2 - begchunk
    if ((ubnd_out - lbnd_out + 1) /= size(vGCM_recvcols)) then
       call endrun(sub//": size of data_out wrong")
    end if

    ! Post receives
    msg_num = 0
    do index = 1, size(vGCM_recvcols)
       if (vGCM_recvcols(index)%phys_tag >= 0) then
          msg_num = msg_num + 1
          call MPI_irecv(data_out(lbnd_out + msg_num - 1), 1, MPI_REAL8,      &
               vGCM_recvcols(index)%phys_task, vGCM_recvcols(index)%phys_tag,    &
               mpicom, recv_request(msg_num), ierr)
       end if
    end do
    ! Post sends
    msg_num = 0
    do index = 1, size(phys_sendcols)
       if (phys_sendcols(index)%phys_tag >= 0) then
          chnk = phys_sendcols(index)%phys_chunk
          col = phys_sendcols(index)%phys_col
          msg_num = msg_num + 1
          call MPI_isend(data_in(col+cooff, chnk+choff), 1, MPI_REAL8,        &
               phys_sendcols(index)%vgcm_task,                                &
               phys_sendcols(index)%phys_tag,                                 &
               mpicom, send_request(msg_num), ierr)
       end if
    end do

    ! Process received messages
    do index = 1, num_vgcm_recvmsgs
       call MPI_waitall(num_vgcm_recvmsgs, recv_request, recv_status, ierr)
    end do

  end subroutine comm_to_vgcm_r8

  subroutine comm_to_vgcm_state(data_in, data_out)
    use cam_abortutils,  only: endrun
    use spmd_utils,      only: mpicom, MPI_STATUS_SIZE

    type(phys_to_vgcm_state_t), intent(in)  :: data_in(:,:)
    type(phys_to_vgcm_state_t), intent(out) :: data_out(:)

    integer :: lbnd_in1, ubnd_in1, lbnd_in2, ubnd_in2, lbnd_out, ubnd_out
    integer :: col, chnk, index
    integer :: choff, cooff ! Offsets in case lower bounds do not match
    integer :: msg_num
    integer :: ierr
    integer :: send_request(num_phys_sendmsgs)
    integer :: recv_request(num_vgcm_recvmsgs)
    integer :: recv_status(MPI_STATUS_SIZE, num_vgcm_recvmsgs)
    character(len=*), parameter :: sub = "comm_to_vgcm_state"

    ! Make sure the dummy arrays match expectations
    lbnd_in1 = LBOUND(data_in, 1)
    ubnd_in1 = UBOUND(data_in, 1)
    lbnd_in2 = LBOUND(data_in, 2)
    ubnd_in2 = UBOUND(data_in, 2)
    lbnd_out = LBOUND(data_out, 1)
    ubnd_out = UBOUND(data_out, 1)
    if ((ubnd_in1 - lbnd_in1 + 1) /= pcols) then
       call endrun(sub//": first index of data_in wrong")
    end if
    cooff = lbnd_in1 - 1
    if ((ubnd_in2 - lbnd_in2) /= (ghost_chunk_e - begchunk)) then
       call endrun(sub//": second index of data_in wrong")
    end if
    choff = lbnd_in2 - begchunk
    if ((ubnd_out - lbnd_out + 1) /= size(vGCM_recvcols)) then
       call endrun(sub//": size of data_out wrong")
    end if

    ! Post receives
    msg_num = 0
    do index = 1, size(vGCM_recvcols)
       if (vGCM_recvcols(index)%phys_tag >= 0) then
          msg_num = msg_num + 1
          call MPI_irecv(data_out(lbnd_out + msg_num - 1), 1,                 &
               vGCM_state_mpi_type, vGCM_recvcols(index)%phys_task,           &
               vGCM_recvcols(index)%phys_tag, mpicom, recv_request(msg_num), ierr)
       end if
    end do
    ! Post sends
    msg_num = 0
    do index = 1, size(phys_sendcols)
       if (phys_sendcols(index)%phys_tag >= 0) then
          chnk = phys_sendcols(index)%phys_chunk
          col = phys_sendcols(index)%phys_col
          msg_num = msg_num + 1
          call MPI_isend(data_in(col+cooff, chnk+choff), 1,                   &
               vGCM_state_mpi_type, phys_sendcols(index)%vgcm_task,           &
               phys_sendcols(index)%phys_tag, mpicom, send_request(msg_num), ierr)
       end if
    end do

    ! Wait until all receive messages complete
    do index = 1, num_vgcm_recvmsgs
       call MPI_waitall(num_vgcm_recvmsgs, recv_request, recv_status, ierr)
    end do

  end subroutine comm_to_vgcm_state

  subroutine comm_to_phys_tend(data_in, data_out_chan1, data_out_chan2)
    use cam_abortutils,  only: endrun
    use spmd_utils,      only: mpicom, MPI_STATUS_SIZE, iam

    type(vgcm_to_phys_tend_t), intent(in)  :: data_in(:)
    type(vgcm_to_phys_tend_t), intent(out) :: data_out_chan1(:,:)
    type(vgcm_to_phys_tend_t), intent(out) :: data_out_chan2(:,:)

    integer :: lbnd_in, ubnd_in, lbnd_out1, ubnd_out1, lbnd_out2, ubnd_out2
    integer :: col, chnk, index
    integer :: choff, cooff ! Offsets in case lower bounds do not match
    integer :: msg_num
    integer :: ierr
    integer :: send_request(num_vgcm_sendmsgs)
    integer :: recv_request(num_phys_recvmsgs)
    integer :: recv_status(MPI_STATUS_SIZE, num_phys_recvmsgs)
    character(len=*), parameter :: sub = "comm_to_phys_tend"

    ! Make sure the dummy arrays match expectations
    if (ANY(SHAPE(data_out_chan1) /= SHAPE(data_out_chan2))) then
       call endrun(sub//": data_out_chan1 and data_out_chan2 must have the same shape")
    end if
    lbnd_in   = LBOUND(data_in, 1)
    ubnd_in   = UBOUND(data_in, 1)
    lbnd_out1 = LBOUND(data_out_chan1, 1)
    ubnd_out1 = UBOUND(data_out_chan1, 1)
    lbnd_out2 = LBOUND(data_out_chan1, 2)
    ubnd_out2 = UBOUND(data_out_chan1, 2)
    if ((ubnd_in - lbnd_in + 1) /= num_vgcm_sendmsgs) then
       call endrun(sub//": size of data_in wrong")
    end if
    if ((ubnd_out1 - lbnd_out1 + 1) /= pcols) then
       call endrun(sub//": first index of data_out wrong")
    end if
    cooff = lbnd_out1 - 1
    if ((ubnd_out2 - lbnd_out2) /= (endchunk - begchunk)) then
       call endrun(sub//": second index of data_out wrong")
    end if
    choff = lbnd_out2 - begchunk

    ! Post receives
    msg_num = 0
    do index = 1, size(phys_sendcols)
       if (phys_sendcols(index)%vgcm_tag >= 0) then
          chnk = phys_sendcols(index)%phys_chunk
          col = phys_sendcols(index)%phys_col
          msg_num = msg_num + 1
          if (phys_sendcols(index)%phys_chan_type == 1) then
             call MPI_irecv(data_out_chan1(col+cooff, chnk+choff), 1,         &
                  vGCM_tend_mpi_type, phys_sendcols(index)%vgcm_task,         &
                  phys_sendcols(index)%vgcm_tag, mpicom, recv_request(msg_num), ierr)
          else
             call MPI_irecv(data_out_chan2(col+cooff, chnk+choff), 1,         &
                  vGCM_tend_mpi_type, phys_sendcols(index)%vgcm_task,         &
                  phys_sendcols(index)%vgcm_tag, mpicom, recv_request(msg_num), ierr)
          end if
       end if
    end do
    if (msg_num /= num_phys_recvmsgs) then
       call endrun("Wrong number of receives posted")
    end if
    ! Post sends
    msg_num = 0
    do index = 1, size(vGCM_recvcols)
       if (vGCM_recvcols(index)%vgcm_tag >= 0) then
          msg_num = msg_num + 1
          call MPI_isend(data_in(lbnd_in + msg_num - 1), 1,                   &
               vGCM_tend_mpi_type, vGCM_recvcols(index)%phys_task,            &
               vGCM_recvcols(index)%vgcm_tag, mpicom, send_request(msg_num), ierr)
       end if
    end do

    ! Wait until all receive messages complete
    do index = 1, num_phys_recvmsgs
       call MPI_waitall(num_phys_recvmsgs, recv_request, recv_status, ierr)
    end do

  end subroutine comm_to_phys_tend

  integer function cons_phys_tag(phys_pe, phys_chunk, phys_col, chan_pe, chan_num)
    use shr_kind_mod, only: i8 => shr_kind_i8
    use spmd_utils,   only: npes
    use q3d_runtime,  only: q3d_total_channels

    integer, intent(in) :: phys_pe
    integer, intent(in) :: phys_chunk
    integer, intent(in) :: phys_col
    integer, intent(in) :: chan_pe
    integer, intent(in) :: chan_num

    integer(kind=i8)    :: temp

    temp = phys_pe
    temp = (temp * npes) + chan_pe
    temp = (temp * ghost_chunk_e) + phys_chunk
    temp = (temp * pcols) + phys_col
    temp = (temp * q3d_total_channels) + chan_num

    cons_phys_tag = int(mod(temp, mpi_tag_max))

  end function cons_phys_tag

  integer function cons_vgcm_tag(chan_pe, chan_num, phys_pe, phys_chunk, phys_col)
    use shr_kind_mod, only: i8 => shr_kind_i8
    use spmd_utils,   only: npes
    use q3d_runtime,  only: q3d_total_channels

    integer, intent(in) :: chan_pe
    integer, intent(in) :: chan_num
    integer, intent(in) :: phys_pe
    integer, intent(in) :: phys_chunk
    integer, intent(in) :: phys_col

    integer(kind=i8)    :: temp

    temp = chan_pe
    temp = (temp * npes) + phys_pe
    temp = (temp * q3d_total_channels) + chan_num
    temp = (temp * pcols) + phys_col
    temp = (temp * ghost_chunk_e) + phys_chunk

    cons_vgcm_tag = int(mod(temp, mpi_tag_max))

  end function cons_vgcm_tag

  integer function q3d_channel_number_to_pe(chan_num) result(pe_num)
    !! Distribute the VVM channels across the available tasks
    !! Clumping adjacent channels on the same task may reduce communication
    !!   requirements (as they are more likely to share ghost columns)
    use spmd_utils,  only: npes
    use q3d_runtime, only: q3d_total_channels

    ! Dummy arg
    integer, intent(in) :: chan_num
    ! Local vars
    integer :: min_chan_num ! Number of channels in every PE
    integer :: plus_pes     ! # PEs with one extra channel
    integer :: overload     ! num of channels in overloaded tasks

    ! Note that chan_num is one based while pes are zero based
    min_chan_num = (q3d_total_channels / npes)
    plus_pes = MOD(q3d_total_channels, npes)
    pe_num = -1
    overload = plus_pes * (min_chan_num + 1)
    if (chan_num > overload) then
       ! This channel will go to a task that is past the overloaded ones
       pe_num = ((chan_num - 1 - overload) / min_chan_num) + plus_pes
    else
       ! This channel will go in an overloaded task
       pe_num = (chan_num - 1) / (min_chan_num + 1)
    end if

  end function q3d_channel_number_to_pe

  subroutine q3d_comm_read_grid(filename)
    use shr_kind_mod,     only: r8=>shr_kind_r8
    use physconst,        only: pi
    use spmd_utils,       only: masterproc, iam, npes, mpicom, MPI_MAX
    use spmd_utils,       only: MPI_INTEGER, MPI_INTEGER8, MPI_ADDRESS_KIND
    use phys_grid,        only: get_ncols_p
    use ppgrid,           only: pcols, begchunk, endchunk
    use cam_grid_support, only: cam_grid_id, cam_grid_dimensions
    use cam_grid_support, only: cam_grid_get_gcid
    use cam_pio_utils,    only: cam_pio_openfile, cam_pio_closefile
    use cam_pio_utils,    only: pio_subsystem
    use pio,              only: file_desc_t, var_desc_t, io_desc_t
    use pio,              only: pio_double, pio_int, PIO_NOERR, PIO_NOWRITE
    use pio,              only: pio_seterrorhandling, PIO_BCAST_ERROR
    use pio,              only: pio_inq_varid, PIO_global
    use pio,              only: pio_inq_dimid, pio_inq_dimlen, PIO_get_att
    use pio,              only: pio_initdecomp, pio_freedecomp, pio_read_darray
    use cam_abortutils,   only: endrun
    use phys_grid,        only: get_rlat_all_p, get_rlon_all_p
    use q3d_runtime,      only: q3d_set_channel_decomp
    use q3d_coupling,     only: q3d_set_indexing

    character(len=*), intent(in)     :: filename

    logical                          :: filefound
    type(file_desc_t)                :: fh
    integer                          :: ncols
    integer                          :: ierr
    integer                          :: err_handling
    ! vGCM grid data
    ! channel vars are allocated pcols*nchunks,halo_len,num chans (2)
    integer,  allocatable            :: source_cols(:)
    integer,  allocatable            :: channel_number(:,:,:)
    integer,  allocatable            :: channel_index(:,:)
    real(r8), allocatable            :: lat(:) ! For lat/lon checks
    real(r8), allocatable            :: lon(:) ! For lat/lon checks
    real(r8)                         :: rlat(pcols) ! For lat/lon checks
    real(r8)                         :: rlon(pcols) ! For lat/lon checks
    ! I/O variables
    integer                          :: dimid
    integer                          :: grid_len  ! grid size dim
    integer                          :: col_len   ! total num columns dim
    integer                          :: col_types ! # col types dim
    integer                          :: halo_len  ! # halo dim
    integer                          :: physgrid
    type(io_desc_t)                  :: iodesc
    type(var_desc_t)                 :: vardesc
    integer                          :: gdims(3)
    integer                          :: i, j, k, h1, h2, gind, lind, crmpe
    integer                          :: coff         ! Channel index offset
    integer                          :: grid_nedge   ! # columns across a face
    integer                          :: cchunk, ccol ! For assigning gcids
    integer                          :: num_chunks
    integer                          :: local_cols ! Total number of local cols
    integer                          :: total_chans, begchan, endchan
    integer(iMap),       pointer     :: gcid(:) ! Global column number
    integer(iMap),       allocatable :: ldof(:) ! Local decomposition
    ! Send/receive variables
    integer                          :: send_cols(0:npes-1) ! # cols to send
    integer                          :: recv_cols(0:npes-1) ! # cols to recv
    integer                          :: send_check(2 * ((nhalo_vGCM * 2) + 1))
    ! MPI variables
    integer                          :: sdispls(0:npes-1)   ! vgcm send offsets
    integer                          :: rdispls(0:npes-1)   ! vgcm recv offsets
    integer                          :: tdispls(0:npes-1)   ! Temp
    integer(kind=MPI_ADDRESS_KIND)   :: offsets(15)   ! For new MPI types
    integer                          :: origtypes(15) ! For new MPI types
    integer                          :: lengths(15)   ! For new MPI types
    integer(kind=MPI_ADDRESS_KIND)   :: extent        ! For new MPI types
    type(phys_to_vgm_loc_t)          :: dummy_loc(2)  ! For new MPI types
    integer                          :: column_type0  ! dummy type
    integer                          :: column_type   ! New MPI type
    ! Convert to degrees to radians
    real(r8),            parameter   :: deg2rad = pi / 180.0_r8
    character(len=128)               :: errmsg
    character(len=*),    parameter   :: sub = 'q3d_comm_read_grid'

    if (masterproc) then
       write(iulog, *) sub, ': Constructing vGCM stencil from ', trim(filename)
       call shr_sys_flush(iulog)
    end if

    ! Retrieve the global column indices for physics columns
    nullify(gcid)
    physgrid = cam_grid_id('physgrid')
    call cam_grid_get_gcid(physgrid, gcid)

    inquire(file=trim(filename), exist=filefound)
    if (.not.filefound) then
       call endrun(sub//': ERROR: could not find grid file '//trim(filename))
    end if
    call cam_pio_openfile(fh, trim(filename), PIO_NOWRITE)
    call pio_seterrorhandling(fh, PIO_BCAST_ERROR, err_handling)

    ierr = PIO_get_att(fh, PIO_global, 'grid_nedge', nedge)

    call cam_grid_dimensions(physgrid, gdims(1:2))

    ! Read in file dimensions
    ierr = pio_inq_dimid(fh, 'grid_size', dimid)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find dimension "grid_size" on '//trim(filename))
    end if
    ierr = pio_inq_dimlen(fh, dimid, grid_len)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR finding length of dimension "grid_size" on '//trim(filename))
    end if
    ! Do some basic grid len checks
    if (grid_len /= PRODUCT(gdims(1:2))) then
       if (masterproc) then
          write(iulog, *) sub, ' ERROR: grid_size = ',grid_len,', expected ',PRODUCT(gdims(1:2))
          call shr_sys_flush(iulog)
       end if
       call endrun(sub//': ERROR: grid size is wrong on '//trim(filename))
    end if
    call MPI_allreduce(MAXVAL(gcid), max_prog_column, 1, MPI_INTEGER8,        &
         MPI_MAX, mpicom, ierr)
    if (max_prog_column /= grid_len) then
       if (masterproc) then
          write(iulog, *) sub, ' ERROR: grid_size = ',grid_len,', max gcid = ',max_prog_column
          call shr_sys_flush(iulog)
       end if
       call endrun(sub//': ERROR: grid size does not match max gcid on '//trim(filename))
    end if
    ! Column size is SE columns plus ghost (off-face) columns
    ierr = pio_inq_dimid(fh, 'column_size', dimid)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find dimension "column_size" on '//trim(filename))
    end if
    ierr = pio_inq_dimlen(fh, dimid, col_len)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR finding length of dimension "column_size" on '//trim(filename))
    end if
    if (col_len < grid_len) then
       write(errmsg, '(a,i0)') 'ERROR: column_size must be at least as large as grid_size'
       call endrun(sub//trim(errmsg))
    end if
    if ((col_len - grid_len) /= (6 * ((nedge * 8) + 12))) then
       write(errmsg, '(2(a,i0))') 'Incorrect number of ghost columns (',      &
            (col_len - grid_len),' should be ',(6 * ((nedge * 8) + 12))
       call endrun(sub//trim(errmsg))
    end if
    ! There are two bands of columns across each vGCM cell (physics column)
    ierr = pio_inq_dimid(fh, 'column_channels', dimid)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find dimension "column_channels" on '//trim(filename))
    end if
    ierr = pio_inq_dimlen(fh, dimid, col_types)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR finding length of dimension "column_channels" on '//trim(filename))
    end if
    if (col_types /= 2) then
       write(errmsg, '(a,i0)') 'ERROR: column_channels must be 2, found ',col_types
       call endrun(sub//trim(errmsg))
    end if
    ! halo_size includes the prognostic column plus the halo on both sides
    ierr = pio_inq_dimid(fh, 'halo_size', dimid)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find dimension "halo_len" on '//trim(filename))
    end if
    ierr = pio_inq_dimlen(fh, dimid, halo_len)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR finding length of dimension "halo_len" on '//trim(filename))
    end if
    ! Make sure halo_len is sensible
    if (halo_len /= (2 * nhalo_vGCM) + 1) then
       write(errmsg, '(a,i0)') 'ERROR: halo_len must be (2 * nhalo_vGCM) + 1, found ',halo_len
       call endrun(sub//trim(errmsg))
    end if

    ! Compute decomp for number of ghost columns
    local_cols = 0 ! Total on-face columns on this task
    num_chunks = endchunk - begchunk + 1
    allocate(ldof(num_chunks * pcols))
    gind = 0
    lind = 0
    do j = 1, num_chunks
       do i = 1, pcols
          gind = gind + 1 ! gcid index
          lind = lind + 1
          if (gcid(gind) > 0) then
             ldof(lind) = gcid(gind)
             local_cols = local_cols + 1
          else
             ldof(lind) = 0
          end if
       end do
    end do

    ! Allocate an array to be used to fill in on-face cells during d_p_coupling
    if (allocated(element_fill_index)) then
       deallocate(element_fill_index)
    end if
    allocate(element_fill_index(local_cols))
    do i = 1, local_cols
       element_fill_index%gcid = 0
       element_fill_index%ie   = 0
       element_fill_index%ioff = 0
       element_fill_index%joff = 0
    end do

    ! Prep for number of ghost columns associated with each grid column
    call PIO_InitDecomp(pio_subsystem, PIO_INT, (/grid_len/), &
         ldof, iodesc)
    ierr = pio_inq_varid(fh, 'num_source_cols', vardesc)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find num_source_cols on '//trim(filename))
    end if
    num_chunks = (endchunk - begchunk + 1)
    allocate(source_cols(pcols * num_chunks))
    source_cols = 0
    call PIO_read_darray(fh, vardesc, iodesc, source_cols, ierr)
    ! Cleanup
    deallocate(ldof)
    call PIO_FreeDecomp(pio_subsystem, iodesc)

    local_cols = SUM(source_cols) ! Total number of ghost columns on this task
    ! Allocate an array used to fill ghost cells during d_p_coupling
    if (allocated(ghost_fill_index)) then
       deallocate(ghost_fill_index)
    end if
    allocate(ghost_fill_index(local_cols))
    do i = 1, local_cols
       ghost_fill_index(i)%gcid = 0
       ghost_fill_index(i)%ghost_id = 0
       ghost_fill_index(i)%ie = 0
       ghost_fill_index(i)%ioff = 0
       ghost_fill_index(i)%joff = 0
       ghost_fill_index(i)%ighost = 0
       ghost_fill_index(i)%jghost = 0
       ghost_fill_index(i)%send_c1_indices(:) = 0
       ghost_fill_index(i)%send_c2_indices(:) = 0
    end do

    ! Figure out how many extra chunks we need for the ghost columns
    deallocate(source_cols)
    num_chunks = local_cols / pcols
    if ((num_chunks * pcols) < local_cols) then
       num_chunks = num_chunks + 1
    end if
    ghost_chunk_s = endchunk + 1
    ghost_chunk_e = endchunk + num_chunks
    allocate(ncols_p(begchunk:ghost_chunk_e))
    do j = begchunk, endchunk
       ncols_p(j) = get_ncols_p(j)
    end do
    do j = ghost_chunk_s, ghost_chunk_e
       if (local_cols <= 0) then
          write(errmsg, '(a,i0,a)') ': Ran out of ghost columns with ', &
               (ghost_chunk_e - j + 1), 'chunks to go'
          call endrun(sub//errmsg)
       end if
       ncols_p(j) = MIN(pcols, local_cols)
       local_cols = local_cols - ncols_p(j)
    end do
    num_chunks = ghost_chunk_e - begchunk + 1
    local_cols = pcols * num_chunks

    ! Allocate column info and set local info
    if (allocated(phys_cols)) then
       deallocate(phys_cols)
    end if
    allocate(phys_cols(pcols, begchunk:ghost_chunk_e))
    gind = 0
    do j = begchunk, endchunk
       do i = 1, pcols
          gind = gind + 1 ! gcid index
          phys_cols(i, j)%gcid = gcid(gind)
          if (gcid(gind) > 0) then
             phys_cols(i, j)%task = iam
          else
             phys_cols(i, j)%task = -1
          end if
       end do
    end do

    ! Now we need to assign all of the ghost chunks, however, we do
    ! not yet know the correct gcids so we have to read them in.
    ! Compute map for ghost source columns (2-D halo)
    allocate(ldof((endchunk-begchunk+1) * pcols * halo_len * halo_len))
    lind = 0
    do h2 = 1, halo_len
       do h1 = 1, halo_len
          gind = 0
          do j = begchunk, endchunk
             do i = 1, pcols
                gind = gind + 1 ! gcid index
                lind = lind + 1
                if (gcid(gind) > 0) then
                   ldof(lind) = ((((h2 - 1) * halo_len) + (h1 - 1)) * grid_len) + gcid(gind)
                else
                   ldof(lind) = 0
                end if
             end do
          end do
       end do
    end do

    call PIO_InitDecomp(pio_subsystem, PIO_INT, (/grid_len, halo_len, halo_len/), ldof, iodesc)
    ! Read in grid_channel_number
    ierr = pio_inq_varid(fh, 'source_columns', vardesc)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find source_columns on '//trim(filename))
    end if
    allocate(channel_number((endchunk-begchunk+1) * pcols, halo_len, halo_len))
    channel_number = -1
    call PIO_read_darray(fh, vardesc, iodesc, channel_number, ierr)
    ! Now, assign all the found ghost columns a chunk and column
    ! Also store info in host column (dynamics source column for ghost info)
    cchunk = endchunk
    ccol = pcols
    do j = begchunk, endchunk
       ncols = ncols_p(j)
       do i = 1, ncols
          do h2 = 1, halo_len
             do h1 = 1, halo_len
                k = ((j - begchunk) * pcols) + i
                phys_cols(i, j)%ghost_source(h1-nhalo_vGCM-1,h2-nhalo_vGCM-1) = channel_number(k, h1, h2)
                if (channel_number(k, h1, h2) > 0) then
                   ccol = ccol + 1
                   if (ccol > ncols_p(cchunk)) then
                      ccol = 1
                      cchunk = cchunk + 1
                      if (cchunk > ghost_chunk_e) then
                         write(errmsg, '(a,i0)') ': Exceeded ghost chunks on ',iam
                         call endrun(sub//trim(errmsg))
                      end if
                   end if
                   phys_cols(ccol, cchunk)%gcid = channel_number(k,h1,h2)
                   phys_cols(ccol, cchunk)%task = iam
                end if
             end do
          end do
       end do
    end do
    ! Check our math
    if (cchunk /= ghost_chunk_e) then
       write(errmsg, '(2(a,i0))') ': Failed to fill ',(ghost_chunk_e - cchunk),&
            ' ghost chunks on ',iam
       call endrun(sub//trim(errmsg))
    end if
    ! Cleanup
    deallocate(channel_number)
    deallocate(ldof)
    call PIO_FreeDecomp(pio_subsystem, iodesc)
    ! Stub out any unused ghost columns
    do j = ghost_chunk_s, ghost_chunk_e
       ncols = ncols_p(j)
       do i = ncols + 1, pcols
          phys_cols(i, j)%gcid = 0
          phys_cols(i, j)%task = -1
       end do
    end do
    ! Now, we need to redo gcid to include all the ghost columns
    deallocate(gcid)
    allocate(gcid((ghost_chunk_e - begchunk + 1) * pcols))
    gind = 0
    do j = begchunk, ghost_chunk_e
       do i = 1, pcols
          gind = gind + 1
          gcid(gind) = phys_cols(i, j)%gcid
       end do
    end do

    ! Prep for lat, lon, and face numbers
    allocate(ldof(local_cols))
    lind = 0
    do j = begchunk, ghost_chunk_e
       do i = 1, pcols
          gind = phys_cols(i, j)%gcid
          lind = lind + 1
          if (gind > 0) then
             ldof(lind) = gind
          else
             ldof(lind) = 0
          end if
       end do
    end do
    call PIO_InitDecomp(pio_subsystem, PIO_DOUBLE, (/col_len/), ldof, iodesc)
    allocate(lat(local_cols))
    allocate(lon(local_cols))

    ! Read in grid_center_lon
    ierr = pio_inq_varid(fh, 'column_lon', vardesc)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find column_lon on '//trim(filename))
    end if
    call PIO_read_darray(fh, vardesc, iodesc, lon, ierr)
    ! Read in grid_center_lat
    ierr = pio_inq_varid(fh, 'column_lat', vardesc)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find column_lat on '//trim(filename))
    end if
    call PIO_read_darray(fh, vardesc, iodesc, lat, ierr)
    ! While re read all the coordinates, we can only check the grid coords
    do j = begchunk, endchunk
       ncols = get_ncols_p(j)
       call get_rlat_all_p(j, pcols, rlat)
       call get_rlon_all_p(j, pcols, rlon)
       do i = 1, ncols
          k = ((j - begchunk) * pcols) + i
          ! Check coordinates against physics
          if (abs(rlat(i) - (lat(k) * deg2rad)) > 1.0e-8) then
             write(errmsg, '(2(a,f7.2))') ': grid center lat, ',lat(k)*deg2rad,' does not match phys lat, ',rlat(i)
             call endrun(sub//trim(errmsg))
          end if
          if (abs(rlon(i) - (lon(k) * deg2rad)) > 1.0e-8) then
             write(errmsg, '(2(a,f7.2))') ': grid center lon, ',lon(k)*deg2rad,' does not match phys lon, ',rlon(i)
             call endrun(sub//trim(errmsg))
          end if
          ! Coords are okay, store them
          phys_cols(i, j)%rlat = lat(k) * deg2rad
          phys_cols(i, j)%rlon = lon(k) * deg2rad
       end do
    end do
    ! Assign the ghost column coordinates
    do j = ghost_chunk_s, ghost_chunk_e
       ncols = ncols_p(j)
       do i = 1, ncols
          phys_cols(i, j)%rlat = lat(k) * deg2rad
          phys_cols(i, j)%rlon = lon(k) * deg2rad
       end do
    end do
    ! Cleanup
    deallocate(lat)
    deallocate(lon)
    call PIO_FreeDecomp(pio_subsystem, iodesc)

    ! Read in the column faces using the same decomposition
    call PIO_InitDecomp(pio_subsystem, PIO_INT, (/col_len/), ldof, iodesc)

    ! Read in column face numbers
    ierr = pio_inq_varid(fh, 'column_face', vardesc)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find column_face on '//trim(filename))
    end if
    allocate(source_cols(local_cols))
    call PIO_read_darray(fh, vardesc, iodesc, source_cols, ierr)
    ! Assign the column face number
    do j = begchunk, ghost_chunk_e
       ncols = ncols_p(j)
       do i = 1, ncols
          phys_cols(i, j)%face = source_cols(((j - begchunk)*pcols) + i)
       end do
    end do
    ! Cleanup
    deallocate(source_cols)
    deallocate(ldof)
    call PIO_FreeDecomp(pio_subsystem, iodesc)

    ! Compute map for column index (no halo)
    allocate(ldof(local_cols * col_types))
    lind = 0
    do k = 1, col_types
       gind = 0
       do j = begchunk, ghost_chunk_e
          do i = 1, pcols
             gind = gind + 1 ! gcid index
             lind = lind + 1
             if (gcid(gind) > 0) then
                ldof(lind) = ((k - 1) * col_len) + gcid(gind)
             else
                ldof(lind) = 0
             end if
          end do
       end do
    end do

    call PIO_InitDecomp(pio_subsystem, PIO_INT, (/col_len, col_types/), &
         ldof, iodesc)
    ! Read in grid_channel_index
    ierr = pio_inq_varid(fh, 'grid_channel_index', vardesc)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find grid_channel_index on '//trim(filename))
    end if
    allocate(channel_index(local_cols, col_types))
    channel_index = -1
    call PIO_read_darray(fh, vardesc, iodesc, channel_index, ierr)
    do j = begchunk, ghost_chunk_e
       ncols = ncols_p(j)
       do i = 1, ncols
          k = ((j - begchunk) * pcols) + i
          phys_cols(i, j)%channel1_index = channel_index(k, 1)
          phys_cols(i, j)%channel2_index = channel_index(k, 2)
       end do
    end do
    ! Cleanup
    deallocate(channel_index)
    deallocate(ldof)
    call PIO_FreeDecomp(pio_subsystem, iodesc)

    ! Compute map for column number (1-D halo)
    allocate(ldof(local_cols * halo_len * col_types))
    lind = 0
    do k = 1, col_types
       do h1 = 1, halo_len
          gind = 0
          do j = begchunk, ghost_chunk_e
             do i = 1, pcols
                gind = gind + 1 ! gcid index
                lind = lind + 1
                if (gcid(gind) > 0) then
                   ldof(lind) = ((((k - 1) * halo_len) + (h1 - 1)) * col_len) + gcid(gind)
                else
                   ldof(lind) = 0
                end if
             end do
          end do
       end do
    end do

    call PIO_InitDecomp(pio_subsystem, PIO_INT, (/col_len, halo_len, col_types/), ldof, iodesc)
    ! Read in grid_channel_number
    ierr = pio_inq_varid(fh, 'grid_channel_number', vardesc)
    if (ierr /= PIO_NOERR) then
       call endrun(sub//': ERROR: could not find grid_channel_number on '//trim(filename))
    end if
    allocate(channel_number(local_cols, halo_len, col_types))
    channel_number = -1
    call PIO_read_darray(fh, vardesc, iodesc, channel_number, ierr)
    do j = begchunk, ghost_chunk_e
       ncols = ncols_p(j)
       do i = 1, ncols
          do h1 = 1, halo_len
             k = ((j - begchunk) * pcols) + i
             phys_cols(i, j)%channel1_num(h1-nhalo_vGCM-1) = channel_number(k,h1,1)
             phys_cols(i, j)%channel2_num(h1-nhalo_vGCM-1) = channel_number(k,h1,2)
          end do
       end do
    end do
    ! Do some checks
    if ((nedge * 3) < MAXVAL(channel_number)) then
       write(errmsg, '(2(a,i0))') ': nedge error, nedge = ',nedge,' but max channel number = ',MAXVAL(channel_number)
       call endrun(sub//trim(errmsg))
    end if
    ! Cleanup
    deallocate(channel_number)
    deallocate(ldof)
    call PIO_FreeDecomp(pio_subsystem, iodesc)

    ! Back to old error handling
    call pio_seterrorhandling(fh, err_handling)
    ! Clean up
    call cam_pio_closefile(fh)
    ! We no longer need gcid
    deallocate(gcid)
    nullify(gcid)

    ! How many channels are on each task?
    ! Maybe we should include the number of cells as a global attribute
    grid_nedge = FLOOR(SQRT(REAL(grid_len,kind=r8)/6._r8) + 1.0E-8_r8)
    if (grid_nedge * grid_nedge * 6 /= grid_len) then
       write(errmsg, '(2(a,i0))') ' ERROR: Bad size len, ',grid_len,', got edge size of ',grid_nedge
       call endrun(sub//trim(errmsg))
    end if
    total_chans = grid_nedge * 3
    ! Distribute channels in chunks of approximately the same size
    k = total_chans / npes ! Number of channels in every PE
    j = MOD(total_chans, npes) ! # PEs with one extra channel
    begchan = (iam * k) + MIN(iam, j) + 1
    endchan = begchan + k - 1
    if (iam < j) then
       endchan = endchan + 1
    end if
    if (masterproc) then
       write(iulog, '(2a,i0)') sub, ': Total VVM channels = ',total_chans
       call shr_sys_flush(iulog)
    end if
    call q3d_set_channel_decomp(total_chans, begchan, endchan)
    ! Figure out how many columns we will be sending to each PE
    !!XXgoldyXX: To do: Implement separate path for crm pe == iam
    send_cols = 0
    do j = begchunk, ghost_chunk_e
       ncols = ncols_p(j)
       do i = 1, ncols
          do h1 = -nhalo_vGCM, nhalo_vGCM
             gind = phys_cols(i, j)%channel1_num(h1)
             if (gind > 0) then
                lind = q3d_channel_number_to_pe(gind)
                if ((lind < 0) .or. (lind >= npes)) then
                   write(errmsg, '(4(a,i0))') ': pe_num (', lind,             &
                        ') out of range for channel1 = ', gind,               &
                        ' for j = ', j-begchunk, ', i = ', i
                   call endrun(sub//trim(errmsg))
                end if
                send_cols(lind) = send_cols(lind) + 1
             end if
             gind = phys_cols(i, j)%channel2_num(h1)
             if (gind > 0) then
                lind = q3d_channel_number_to_pe(gind)
                if ((lind < 0) .or. (lind >= npes)) then
                   write(errmsg, '(4(a,i0))') ': pe_num (', lind,             &
                        ') out of range for channel2 = ', gind,               &
                        ' for j = ', j-begchunk, ', i = ', i
                   call endrun(sub//trim(errmsg))
                end if
                send_cols(lind) = send_cols(lind) + 1
             end if
          end do
       end do
    end do
    ! Send phys column numbers to vGCM PEs
    call MPI_alltoall(send_cols, 1, MPI_INTEGER, recv_cols, 1, MPI_INTEGER,   &
         mpicom, ierr)

    ! Allocate memory for sending and receiving columns
    allocate(vGCM_recvcols(SUM(recv_cols)))
    allocate(phys_sendcols(SUM(send_cols)))
    do i = 1, size(phys_sendcols)
       phys_sendcols(i)%phys_task = -1
    end do
    num_phys_sendmsgs = 0

    ! Fill in the sending side
    ! First, turn send_cols into displacements
    sdispls(0) = 0
    do i = 1, npes - 1
       sdispls(i) = sdispls(i - 1) + send_cols(i - 1)
    end do
    ! tdispls is used to fill in the send columns
    tdispls(:) = sdispls(:)
    if ((tdispls(npes - 1) + send_cols(npes-1)) /= size(phys_sendcols)) then
       write(errmsg, *) ': ',tdispls(npes - 1),send_cols(npes-1),size(phys_sendcols)
       call endrun(sub//trim(errmsg))
    end if
    ! Now, turn recv_cols into displacements
    rdispls(0) = 0
    do i = 1, npes - 1
       rdispls(i) = rdispls(i - 1) + recv_cols(i - 1)
    end do
    ! Find the ghost id offset for each source col
    lind = 0
    allocate(ghost_index(pcols, ghost_chunk_s:ghost_chunk_e))
    ghost_index = 0
    do j = begchunk, endchunk
       do i = 1, ncols_p(j)
          do h2 = -nhalo_vGCM, nhalo_vGCM
             do h1 = -nhalo_vGCM, nhalo_vGCM
                if (phys_cols(i,j)%ghost_source(h1, h2) > 0) then
                   if ((h1 == 0) .and. (h2 == 0)) then
                      call endrun(sub//': Error: ghost_source on source column')
                   end if
                   lind = lind + 1
                   ghost_fill_index(lind)%gcid = phys_cols(i,j)%gcid
                   ghost_fill_index(lind)%ghost_id = phys_cols(i,j)%ghost_source(h1, h2)
                   ghost_fill_index(lind)%ighost = h1
                   ghost_fill_index(lind)%jghost = h2
                   ! We need a way to find this later
                   filefound = .false.
                   do cchunk = ghost_chunk_s, ghost_chunk_e
                      do ccol = 1, ncols_p(cchunk)
                         if (phys_cols(ccol, cchunk)%gcid == ghost_fill_index(lind)%ghost_id) then
                            if (ghost_index(ccol, cchunk) > 0) then
                               call endrun(sub//': Error: ghost_index slot already used')
                            end if
                            ghost_index(ccol, cchunk) = lind
                            filefound = .true.
                            exit
                         end if
                      end do
                      if (filefound) then
                         exit
                      end if
                   end do
                   if (.not. filefound) then
                      call endrun(sub//': Did not find ghost cell index')
                   end if
                end if
             end do
          end do
       end do
    end do
    if (lind < size(ghost_fill_index)) then
       write(errmsg, '(a,i0,a)') ': ghost_fill_index has ',size(ghost_fill_index)-lind,' unused slots'
       call endrun(sub//trim(errmsg))
    end if
    ! Find the element id offset for each on-face column
    lind = 0
    allocate(element_index(pcols, begchunk:endchunk))
    element_index = 0
    do j = begchunk, endchunk
       do i = 1, ncols_p(j)
          lind = lind + 1
          element_fill_index(lind)%gcid = phys_cols(i,j)%gcid
          element_index(i, j) = lind
       end do
    end do
    ! Now do the fill
    do j = begchunk, ghost_chunk_e
       ncols = ncols_p(j)
       do i = 1, ncols
          k = 1
          send_check = -1
          ! Loop over all the channels where this column could be sent
          ! zero is its prognostic channel
          ! Some slots will be unused (col is within nhalo_vGCM of a face edge)
          do h1 = -nhalo_vGCM, nhalo_vGCM
             ! Channel 1
             gind = phys_cols(i, j)%channel1_num(h1)
             if (gind > 0) then
                crmpe = q3d_channel_number_to_pe(gind)
                tdispls(crmpe) = tdispls(crmpe) + 1
                h2 = tdispls(crmpe)
                if (phys_sendcols(h2)%phys_task /= -1) then
                   write(errmsg, *) ': phys_sendcols(',h2,') = ',phys_sendcols(h2)%phys_task,', crmpe = ',crmpe
                   call endrun(sub//trim(errmsg))
                end if
                phys_sendcols(h2)%phys_task = iam
                phys_sendcols(h2)%phys_chunk = j
                phys_sendcols(h2)%phys_col = i
                phys_sendcols(h2)%phys_chan_type = 1
                phys_sendcols(h2)%vgcm_task = crmpe
                phys_sendcols(h2)%vgcm_channel_num = gind
                phys_sendcols(h2)%vgcm_channel_seg  = get_channel_segment(gind, phys_cols(i, j)%face)
                if (phys_sendcols(h2)%vgcm_channel_seg < 1) then
                   write(errmsg, '(3(a,i0))') ': bad channel1 segment for ',gind,", ",phys_cols(i, j)%channel1_index
                   call endrun(sub//trim(errmsg))
                end if
                ! vgcm_channel_index is relative to the channel segment
                coff = (phys_sendcols(h2)%vgcm_channel_seg - 1) * grid_nedge
                phys_sendcols(h2)%vgcm_channel_index = phys_cols(i, j)%channel1_index - coff
                phys_sendcols(h2)%vgcm_ghost_index = h1
                if (bad_ghost_column(phys_cols(i, j)%face, 1, gind,           &
                     phys_sendcols(h2)%vgcm_channel_index, h1, grid_nedge)) then
                   write(errmsg, '(a,5(i2,", "),i0)') ': Bad ghost1 set: ',   &
                        gind, phys_sendcols(h2)%vgcm_channel_seg, h1,         &
                        phys_cols(i, j)%channel1_index,coff,                  &
                        phys_cols(i, j)%face
                   call endrun(sub//trim(errmsg))
                end if
                ! We only have to send a column to a vGCM task once
                if (ANY(send_check(:) == phys_sendcols(h2)%vgcm_task)) then
                   phys_sendcols(h2)%phys_tag = -1
                else
                   phys_sendcols(h2)%phys_tag = cons_phys_tag(iam, j, i, crmpe, gind)
                   send_check(k) = phys_sendcols(h2)%vgcm_task
                   num_phys_sendmsgs = num_phys_sendmsgs + 1
                   k = k + 1
                end if
                ! We only send prognostic columns back to physics
                if ((h1 == 0) .and. (j <= endchunk)) then
                   phys_sendcols(h2)%vgcm_tag = cons_vgcm_tag(crmpe, gind, iam, j, i)
                else
                   phys_sendcols(h2)%vgcm_tag = -1
                end if
                if (j >= ghost_chunk_s) then
                   ! Set up a ghost_fill_index slot
                   lind = ghost_index(i, j)
                   if (ghost_fill_index(lind)%send_c1_indices(h1) > 0) then
                      write(errmsg, '(5(a,i0))')                              &
                           ': Already have c1 ghost send index ',h1,          &
                           ' filled at ',iam,', ',i,', ',j,', ',lind
                      call endrun(sub//trim(errmsg))
                   end if
                   ghost_fill_index(lind)%send_c1_indices(h1) = h2
                end if
             end if
             ! Channel 2
             gind = phys_cols(i, j)%channel2_num(h1)
             if (gind > 0) then
                crmpe = q3d_channel_number_to_pe(gind)
                tdispls(crmpe) = tdispls(crmpe) + 1
                h2 = tdispls(crmpe)
                phys_sendcols(h2)%phys_task = iam
                phys_sendcols(h2)%phys_chunk = j
                phys_sendcols(h2)%phys_col = i
                phys_sendcols(h2)%phys_chan_type = 2
                phys_sendcols(h2)%vgcm_task = crmpe
                phys_sendcols(h2)%vgcm_channel_num = gind
                phys_sendcols(h2)%vgcm_channel_seg  = get_channel_segment(gind, phys_cols(i, j)%face)
                if (phys_sendcols(h2)%vgcm_channel_seg < 1) then
                   write(errmsg, '(3(a,i0))') ': bad channel2 segment for ',gind,", ",phys_cols(i, j)%channel2_index
                   call endrun(sub//trim(errmsg))
                end if
                ! vgcm_channel_index is relative to the channel segment
                coff = (phys_sendcols(h2)%vgcm_channel_seg - 1) * grid_nedge
                phys_sendcols(h2)%vgcm_channel_index = phys_cols(i, j)%channel2_index - coff
                phys_sendcols(h2)%vgcm_ghost_index = h1
                if (bad_ghost_column(phys_cols(i, j)%face, 2, gind,           &
                     phys_sendcols(h2)%vgcm_channel_index, h1, grid_nedge)) then
                   write(errmsg, '(a,5(i2,", "),i0)') ': Bad ghost2 set: ',   &
                        gind, phys_sendcols(h2)%vgcm_channel_seg, h1,         &
                        phys_cols(i, j)%channel2_index,coff,                  &
                        phys_cols(i, j)%face
                   call endrun(sub//trim(errmsg))
                end if
                ! We only have to send a column to a vGCM task once
                if (ANY(send_check(:) == phys_sendcols(h2)%vgcm_task)) then
                   phys_sendcols(h2)%phys_tag = -1
                else
                   phys_sendcols(h2)%phys_tag = cons_phys_tag(iam, j, i, crmpe, gind)
                   send_check(k) = phys_sendcols(h2)%vgcm_task
                   num_phys_sendmsgs = num_phys_sendmsgs + 1
                   k = k + 1
                end if
                ! We only send prognostic columns back to physics
                if ((h1 == 0) .and. (j <= endchunk)) then
                   phys_sendcols(h2)%vgcm_tag = cons_vgcm_tag(crmpe, gind, iam, j, i)
                else
                   phys_sendcols(h2)%vgcm_tag = -1
                end if
                if (j >= ghost_chunk_s) then
                   ! Set up a ghost_fill_index slot
                   lind = ghost_index(i, j)
                   if (ghost_fill_index(lind)%send_c2_indices(h1) > 0) then
                      write(errmsg, '(5(a,i0))')                              &
                           ': Already have c2 ghost send index ',h1,          &
                           ' filled at ',iam,', ',i,', ',j,', ',lind
                      call endrun(sub//trim(errmsg))
                   end if
                   ghost_fill_index(lind)%send_c2_indices(h1) = h2
                end if
             end if
          end do
       end do
    end do
    ! We should have used up every slot
    do i = 0, npes - 2
       if (tdispls(i) < sdispls(i + 1)) then
          write(errmsg, '(a,i0,a)') ': tdispls(',i,') slot not filled'
          call endrun(sub//':  Not all tdispls slots filled')
       else if (tdispls(i) > sdispls(i + 1)) then
          write(errmsg, '(a,i0,a)') ': tdispls(',i,') slot over filled'
          call endrun(sub//':  Not all tdispls slots filled')
       end if
    end do
    if (tdispls(npes - 1) < size(phys_sendcols)) then
       write(errmsg, '(a,i0,a)') ':  tdispls(',i,') slot not filled'
       call endrun(sub//trim(errmsg))
    else if (tdispls(npes - 1) > size(phys_sendcols)) then
       write(errmsg, '(a,i0,a)') ':  tdispls(',i,') slot over filled'
       call endrun(sub//trim(errmsg))
    end if
    ! Fill in the rest of ghost_fill_index
    call q3d_set_indexing(ghost_fill_index)
    ! Fill in the rest of element_fill_index
    call q3d_set_indexing(element_fill_index)

    ! Create MPI data type for columns
    lengths(:) = 1
    origtypes(:) = MPI_INTEGER
    h1 = 0
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%phys_task,          offsets(h1), ierr)
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%phys_chunk,         offsets(h1), ierr)
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%phys_col,           offsets(h1), ierr)
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%phys_chan_type,     offsets(h1), ierr)
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%vgcm_task,          offsets(h1), ierr)
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%vgcm_channel_num,   offsets(h1), ierr)
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%vgcm_channel_seg,   offsets(h1), ierr)
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%vgcm_channel_index, offsets(h1), ierr)
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%vgcm_ghost_index,   offsets(h1), ierr)
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%phys_tag,           offsets(h1), ierr)
    h1 = h1 + 1
    call MPI_Get_address(dummy_loc(1)%vgcm_tag,           offsets(h1), ierr)
    do i = h1, 1, -1
       offsets(i) = offsets(i) - offsets(1)
    end do
    call MPI_type_create_struct(h1, lengths(1:h1), offsets(1:h1), origtypes(1:h1), column_type0, ierr)
    ! Adjust for padding
    call MPI_Get_address(dummy_loc(1)%phys_task, offsets(1), ierr)
    call MPI_Get_address(dummy_loc(2)%phys_task, offsets(2), ierr)
    extent = offsets(2) - offsets(1)
    call MPI_type_create_resized(column_type0, 0_MPI_ADDRESS_KIND, extent, &
         column_type, ierr)
    call MPI_type_commit(column_type, ierr)

    ! Send vGCM to CRM info to PEs
    call MPI_alltoallv(phys_sendcols, send_cols, sdispls, column_type,        &
         vGCM_recvcols, recv_cols, rdispls, column_type, mpicom, ierr)
    if (masterproc) then
       write(iulog, *) sub, ": Finished sending column to CRM info blocks"
       call shr_sys_flush(iulog)
    end if

    ! Fill in the receiving side
    num_vgcm_recvmsgs = 0
    num_vgcm_sendmsgs = 0
    send_cols = 0
    allocate(recv_indexing(size(vGCM_recvcols)))
    i = -1
    j = -1
    if (size(vGCM_recvcols) > 0) then
       if (vGCM_recvcols(1)%phys_tag < 0) then
          call endrun(sub//'Bad start to vGCM_recvcols')
       end if
    end if
    do h1 = 1, size(vGCM_recvcols)
       if (vGCM_recvcols(h1)%phys_tag >= 0) then
          num_vgcm_recvmsgs = num_vgcm_recvmsgs + 1
          gind = vGCM_recvcols(h1)%phys_task
          j = vGCM_recvcols(h1)%phys_chunk
          i = vGCM_recvcols(h1)%phys_col
          h2 = h1
       else if ((vGCM_recvcols(h1)%phys_task /= gind) .or. &
            (vGCM_recvcols(h1)%phys_chunk /= j)       .or. &
            (vGCM_recvcols(h1)%phys_col /= i)) then
          ! The last recvcol does not match this one
          write(errmsg, '(2(a,i0))') ':  No recv index for ',iam,', ',h1
          call endrun(sub//trim(errmsg))
          ! else everything is okay, we just have a repeat to limit MPI messages
       end if
       recv_indexing(h1) = h2
       if (vGCM_recvcols(h1)%vgcm_tag >= 0) then
          num_vgcm_sendmsgs = num_vgcm_sendmsgs + 1
          send_cols(gind) = send_cols(gind) + 1
       end if
    end do
    ! Send vGCM prognostic column numbers to phys PEs
    recv_cols = 0
    call MPI_alltoall(send_cols, 1, MPI_INTEGER, recv_cols, 1, MPI_INTEGER,   &
         mpicom, ierr)
    num_phys_recvmsgs = SUM(recv_cols)

  end subroutine q3d_comm_read_grid

  integer function get_channel_segment(channel_number, face)
    use cam_abortutils,  only: endrun

    integer, intent(in) :: channel_number
    integer, intent(in) :: face

    integer :: channel_type
    integer :: channel_set

    if ((face < 1) .or. (face > 6)) then
       call endrun('get_channel_segment: Bad face number')
    end if
    channel_set = ((channel_number - 1) / nedge) + 1
    if (channel_set == 1) then
       channel_type = 1
    else if (channel_set == 2) then
       if (face >= 5) then
          channel_type = 1
       else
          channel_type = 2
       end if
    else
       channel_type = 2
    end if
    if (channel_number < 1) then
       get_channel_segment = -1
    else
       get_channel_segment = channel_segment(face, channel_type)
    end if
  end function get_channel_segment

  logical function bad_ghost_column(face, chan_type,                          &
       chan_num, chan_seg_ind, ghost_ind, grid_nedge) result(bad)
    use cam_abortutils,  only: endrun

    integer, intent(in) :: face
    integer, intent(in) :: chan_type
    integer, intent(in) :: chan_num
    integer, intent(in) :: chan_seg_ind
    integer, intent(in) :: ghost_ind
    integer, intent(in) :: grid_nedge

    logical             :: index_on_edge
    integer             :: min_chan_num
    integer             :: max_chan_num

    if ((chan_type < 1) .or. (chan_type > 2)) then
       call endrun('bad_ghost_column: bad chan_type')
    end if
    bad = .false.
    index_on_edge = (chan_seg_ind <= -nhalo_vGCM) .or. (chan_seg_ind >= grid_nedge+nhalo_vGCM)
    select case(face)
    case(1,3)
       if (chan_type == 1) then
          min_chan_num = 1
          max_chan_num = grid_nedge
       else
          min_chan_num = grid_nedge + 1
          max_chan_num = 2 * grid_nedge
       end if
    case(2,4)
       if (chan_type == 1) then
          min_chan_num = 1
          max_chan_num = grid_nedge
       else
          min_chan_num = (2 * grid_nedge) + 1
          max_chan_num =3 * grid_nedge
       end if
    case(5,6)
       if (chan_type == 1) then
          min_chan_num = grid_nedge + 1
          max_chan_num = 2 * grid_nedge
       else
          min_chan_num = (2 * grid_nedge) + 1
          max_chan_num =3 * grid_nedge
       end if
    case default
       call endrun('bad_ghost_column: bad face number')
    end select
    if ( (chan_num == min_chan_num) .and.                                     &
         (ghost_ind == -nhalo_vGCM) .and. index_on_edge) then
       bad = .true.
    end if
    if ( (chan_num == max_chan_num) .and.                                     &
         (ghost_ind == nhalo_vGCM) .and. index_on_edge) then
       bad = .true.
    end if
  end function bad_ghost_column

end module q3d_comm
