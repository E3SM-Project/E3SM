module dyn_grid
  !------------------------------------------------------------------------------------------------- 
  ! 
  ! Purpose: Definition of dynamics computational grid.
  ! 
  !          The grid used by the SE dynamics is called the GLL grid, which is 
  !          made up up of elements that correspond to "blocks". The columns 
  !          within each element are located at the Gauss-Lobatto-Legendre (GLL)
  !          quadrature points. The GLL nodes can be used as the physics columns
  !          as well, but this can lead to grid imprinting. The physics can also
  !          be done on a quasi-equal area FV grid when fv_nphys>0, which will 
  !          define fv_nphys^2 cells per element.
  ! 
  ! Method: Variables are private; interface routines used to extract
  !         information for use in user code. Global column index range
  !         defined using full (unreduced) grid. 
  ! 
  ! Entry points:
  !      get_block_bounds_d       get first and last indices in global 
  !                               block ordering
  !      get_block_gcol_d         get column indices for given block
  !      get_block_gcol_cnt_d     get number of columns in given block
  !      get_block_lvl_cnt_d      get number of vertical levels in column
  !      get_block_levels_d       get vertical levels in column
  !      get_gcol_block_d         get global block indices and local columns 
  !                               index for given global column index
  !      get_gcol_block_cnt_d     get number of blocks containing data
  !                               from a given global column index
  !      get_block_owner_d        get process "owning" given block
  !      get_horiz_grid_d         get horizontal grid coordinates and associated
  !                               information
  !      get_horiz_grid_dim_d     get horizontal dimensions of dynamics grid
  !      dyn_grid_get_pref        get reference pressures for the dynamics grid
  !      dyn_grid_get_elem_coords get coordinates of a specified block element 
  !                               of the dynamics grid
  !      dyn_grid_find_gcols      finds nearest column for given lat/lon
  !      dyn_grid_get_colndx      get element block/column and MPI process indices 
  !                               corresponding to a specified global column index
  !
  ! Author: Jim Edwards and Patrick Worley
  ! 
  !-------------------------------------------------------------------------------------------------
  use element_mod,             only: element_t
  use cam_logfile,             only: iulog
  use cam_abortutils,          only: endrun
  use shr_kind_mod,            only: r8 => shr_kind_r8
  use shr_const_mod,           only: SHR_CONST_PI
  use cam_grid_support,        only: iMap
  use dimensions_mod,          only: nelem, nelemd, nelemdmax, ne, np, npsq
  use spmd_utils,              only: iam, mpi_integer, mpi_real8, mpicom, npes, masterproc
  use parallel_mod,            only: par
  use scamMod,                 only: single_column, scm_multcols

  implicit none
  private
  save

  ! The SE dynamics grid
  integer, parameter, public :: dyn_decomp = 101
  integer, parameter, public :: physgrid_d = 102

  ! FV physics grid resolution (physics on GLL grid if NPG=0)
  integer, parameter, public :: fv_nphys = NPG  

  integer                    :: ngcols_d = 0          ! number of dynamics columns
  logical                    :: gblocks_need_initialized = .true.
  integer, public, parameter :: ptimelevels = 2
  real(r8),        parameter :: rad2deg = 180.0_r8 / SHR_CONST_PI

  type block_global_data
     integer :: UniquePtOffset
     integer :: NumUniqueP
     integer :: LocalID
     integer :: Owner
  end type block_global_data

  type(block_global_data), allocatable :: gblocks(:)
  type(element_t), public, pointer :: elem(:) => null()

  public :: dyn_grid_init, define_cam_grids
  public :: get_block_owner_d, get_gcol_block_d, get_gcol_block_cnt_d
  public :: get_block_gcol_cnt_d, get_horiz_grid_dim_d, get_block_levels_d
  public :: get_block_gcol_d, get_block_bounds_d, get_horiz_grid_d
  public :: get_block_lvl_cnt_d, set_horiz_grid_cnt_d, get_dyn_grid_parm
  public :: get_dyn_grid_parm_real2d, get_dyn_grid_parm_real1d
  public :: dyn_grid_get_pref
  public :: dyn_grid_get_elem_coords
  public :: dyn_grid_find_gcols
  public :: dyn_grid_get_colndx
  public :: physgrid_copy_attributes_d
  public :: fv_physgrid_init, fv_physgrid_final
  public :: compute_global_area
  public :: compute_global_coords

  integer(kind=iMap), pointer :: fdofP_local(:,:) => null()

  real(r8), public, pointer :: w(:) => null()       ! weights

  ! Local lat/lon arrays
  real(r8),      allocatable :: pelat_deg(:)        ! pe-local latitudes (degrees)
  real(r8),      allocatable :: pelon_deg(:)        ! pe-local longitudes (degrees)
  real(r8),          pointer :: pearea(:) => null() ! pe-local areas
  integer(iMap),     pointer :: pemap(:)  => null() ! pe-local map for PIO decomp

!===================================================================================================
contains
!===================================================================================================

  subroutine dyn_grid_init()

  end subroutine dyn_grid_init
  !
  !=================================================================================================
  !
  subroutine get_block_bounds_d(block_first,block_last)
    !---------------------------------------------------------------------------                   
    ! Purpose: Return first and last indices used in global block ordering
    ! 
    ! Author: Jim Edwards
    !------------------------------Arguments------------------------------------
    integer, intent(out) :: block_first  ! first (global) index used for blocks
    integer, intent(out) :: block_last   ! last (global) index used for blocks
    !---------------------------------------------------------------------------
    block_first = 1
    block_last = nelem
    return
  end subroutine get_block_bounds_d
  !
  !=================================================================================================
  !
  subroutine get_block_gcol_d(blockid,size,cdex)
    !---------------------------------------------------------------------------                 
    ! Purpose: Return list of dynamics column indices in given block
    ! 
    ! Author: Jim Edwards
    !---------------------------------------------------------------------------
    implicit none
    !------------------------------Arguments------------------------------------
    integer, intent(in ):: blockid      ! global block id
    integer, intent(in ):: size         ! array size
    integer, intent(out):: cdex(size)   ! global column indices
    !----------------------------Local-Variables--------------------------------
    integer :: ic
    !---------------------------------------------------------------------------
    if (fv_nphys > 0) then
      cdex(1) = (blockid-1)*fv_nphys*fv_nphys + 1
      do ic = 2, size
         cdex(ic) = cdex(1) + ic - 1
      end do
    else
      if(gblocks_need_initialized) call gblocks_init()
      do ic=1,size
         cdex(ic)=gblocks(blockid)%UniquePtOffset+ic-1
      end do
    end if ! fv_nphys > 0

    return
  end subroutine get_block_gcol_d
  !
  !=================================================================================================
  !
  integer function get_block_gcol_cnt_d(blockid)
    !---------------------------------------------------------------------------                   
    ! Purpose: Return number of dynamics columns in indicated block
    ! 
    ! Author: Jim Edwards
    !------------------------------Arguments------------------------------------
    integer, intent(in) :: blockid
    !----------------------------Local-Variables--------------------------------
    integer :: ie
    !---------------------------------------------------------------------------
    if (fv_nphys > 0) then
      get_block_gcol_cnt_d = fv_nphys*fv_nphys
    else
      if(gblocks_need_initialized) call gblocks_init()
      get_block_gcol_cnt_d=gblocks(blockid)%NumUniqueP
    end if ! fv_nphys > 0

    return
  end function get_block_gcol_cnt_d
  !
  !=================================================================================================
  !
  integer function get_block_lvl_cnt_d(blockid,bcid)
    use pmgrid, only: plevp
    !---------------------------------------------------------------------------
    ! Purpose: Return number of levels in indicated column. If column
    !          includes surface fields, then it is defined to also
    !          include level 0.
    !
    ! Author: Patrick Worley
    !---------------------------------------------------------------------------
    implicit none
    !------------------------------Arguments------------------------------------
    integer, intent(in) :: blockid  ! global block id
    integer, intent(in) :: bcid    ! column index within block
    !---------------------------------------------------------------------------
    get_block_lvl_cnt_d = plevp
    return
  end function get_block_lvl_cnt_d
  !
  !=================================================================================================
  !
  subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)
    use pmgrid,         only: plev
    use cam_abortutils, only: endrun
    !---------------------------------------------------------------------------
    ! Purpose: Return level indices in indicated column. If column
    !          includes surface fields, then it is defined to also
    !          include level 0.
    !
    ! Author: Patrick Worley
    !---------------------------------------------------------------------------
    implicit none
    !------------------------------Arguments------------------------------------
    integer, intent(in ) :: blockid         ! global block id
    integer, intent(in ) :: bcid            ! column index within block
    integer, intent(in ) :: lvlsiz          ! dimension of levels array
    integer, intent(out) :: levels(lvlsiz)  ! levels indices for block
    !----------------------------Local-Variables--------------------------------
    integer k                        ! loop index
    !---------------------------------------------------------------------------
    !  latitude slice block
    if (lvlsiz < plev + 1) then
       write(iulog,*)'GET_BLOCK_LEVELS_D: levels array not large enough (', &
            lvlsiz,' < ',plev + 1,' ) '
       call endrun
    else
       do k=0,plev
          levels(k+1) = k
       end do
       do k=plev+2,lvlsiz
          levels(k) = -1
       end do
    end if

    return
  end subroutine get_block_levels_d
  !
  !=================================================================================================
  !
  subroutine get_gcol_block_d(gcol, cnt, blockid, bcid, localblockid)
    use kinds,          only: int_kind
    use cam_abortutils, only: endrun
    !---------------------------------------------------------------------------
    ! Purpose: Return global block index and local column index for 
    !          global column index. Element array is naturally ordered.
    ! 
    ! Author: Jim Edwards
    ! Modified: replaced linear search with best guess and binary search   
    !           (Pat Worley, 7/2/09)
    !---------------------------------------------------------------------------
    implicit none
    !------------------------------Arguments------------------------------------
    integer, intent(in ) :: gcol          ! global column index
    integer, intent(in ) :: cnt           ! size of blockid and bcid arrays
    integer, intent(out) :: blockid(cnt)  ! block index
    integer, intent(out) :: bcid(cnt)     ! column index within block
    integer, intent(out), optional :: localblockid(cnt)
    !----------------------------Local-Variables--------------------------------
    integer :: sb, eb, ie, high, low
    logical :: found
    integer :: nphys_sq
    integer, save :: iedex_save = 1
    !---------------------------------------------------------------------------
    if(gblocks_need_initialized) call gblocks_init()
    if (fv_nphys > 0) then

      ! Note: nphys_sq is to fool the compiler in debug mode when fv_nphys=0
      nphys_sq = fv_nphys*fv_nphys
      blockid(1) = 1 + (gcol-1) / nphys_sq
      bcid(1) = 1 + mod(gcol-1, nphys_sq)

      if (present(localblockid)) then
         localblockid = -1
         if (iam == gblocks(blockid(1))%Owner) then
            if (blockid(1) == elem(iedex_save)%globalid) then
               localblockid = iedex_save
            else
               do ie = 1,nelemd
                  if (blockid(1) == elem(ie)%globalid) then
                     localblockid = ie
                     iedex_save = ie
                     exit
                  end if
               end do
            end if
         end if
      end if

    else ! physics is on GLL grid

      found = .false.
      low = 1
      high = nelem

      ! check whether previous found element is the same here
      if (.not. found) then
        ie = iedex_save
        sb = gblocks(ie)%UniquePtOffset
        if (gcol >= sb) then
          eb = sb+gblocks(ie)%NumUniqueP
          if (gcol < eb) then
            found = .true.
          else
            low = ie
          endif
        else
          high = ie
        endif
      endif

      ! check whether next element  is the one wanted
      if ((.not. found) .and. &
           ((low .eq. iedex_save) .or. (iedex_save .eq. nelem))) then
        ie = iedex_save + 1
        if (ie > nelem) ie = 1

        sb = gblocks(ie)%UniquePtOffset
        if (gcol >= sb) then
          eb = sb+gblocks(ie)%NumUniqueP
          if (gcol < eb) then
            found = .true.
          else
            low = ie
          endif
        else
          high = ie
        endif
      endif

      ! otherwise, use a binary search to find element
      if (.not. found) then
        ! (start with a sanity check)
        ie = low
        sb = gblocks(ie)%UniquePtOffset

        ie = high
        eb = gblocks(ie)%UniquePtOffset + gblocks(ie)%NumUniqueP

        if ((gcol < sb) .or.  (gcol >= eb)) then
          do ie=1,nelemd
            write(iulog,*) __LINE__,ie,elem(ie)%idxP%UniquePtOffset,elem(ie)%idxP%NumUniquePts
          end do
          call endrun()
        end if

        do while (.not. found)
          ie = low + (high-low)/2;

          sb = gblocks(ie)%UniquePtOffset
          if (gcol >= sb) then
            eb = sb+gblocks(ie)%NumUniqueP
            if (gcol < eb) then
              found = .true.
            else
              low = ie+1
            endif
          else
            high = ie-1
          endif
        enddo
      endif

      blockid(1)=ie
      if(present(localblockid)) localblockid(1)=gblocks(ie)%LocalID
      bcid(1)=gcol-sb+1
      iedex_save = ie

    end if ! fv_nphys > 0

    return
  end subroutine get_gcol_block_d
  !
  !=================================================================================================
  !
  integer function get_gcol_block_cnt_d(gcol)
    !---------------------------------------------------------------------------
    ! Purpose: Return number of blocks contain data for the vertical column
    !          with the given global column index
    ! 
    ! Author: Patrick Worley
    !---------------------------------------------------------------------------
    implicit none
    !------------------------------Arguments------------------------------------
    integer, intent(in) :: gcol     ! global column index
    !---------------------------------------------------------------------------
    get_gcol_block_cnt_d = 1
    return
  end function get_gcol_block_cnt_d
  !
  !=================================================================================================
  !
  integer function get_block_owner_d(blockid)
    !---------------------------------------------------------------------------                          
    ! Purpose: Return id of processor that "owns" the indicated block
    ! 
    ! Author: Jim Edwards
    !---------------------------------------------------------------------------
    implicit none
    !------------------------------Arguments------------------------------------
    integer, intent(in) :: blockid  ! global block id
    !---------------------------------------------------------------------------
    if (gblocks_need_initialized) call gblocks_init()
    
    if (gblocks(blockid)%Owner>-1) then
      get_block_owner_d=gblocks(blockid)%Owner
    else
      call endrun('Block owner not assigned in gblocks_init')
    end if

   return
  end function get_block_owner_d
  !
  !=================================================================================================
  !
  subroutine get_horiz_grid_dim_d(hdim1_d,hdim2_d)
    !---------------------------------------------------------------------------
    ! Purpose: Returns declared horizontal dimensions of computational grid.
    !          For non-lon/lat grids, declare grid to be one-dimensional,
    !          i.e., (ngcols_d x 1)
    ! 
    ! Author: Patrick Worley
    !------------------------------Arguments------------------------------------
    integer, intent(out)           :: hdim1_d  ! first horizontal dimension
    integer, intent(out), optional :: hdim2_d  ! second horizontal dimension
    !---------------------------------------------------------------------------
    if (fv_nphys > 0) then
      hdim1_d = fv_nphys*fv_nphys*nelem
    else
      hdim1_d = ngcols_d
    end if ! fv_nphys > 0
    
    if (present(hdim2_d)) hdim2_d = 1

    return
  end subroutine get_horiz_grid_dim_d
  !
  !=================================================================================================
  !
  subroutine set_horiz_grid_cnt_d(NumUniqueCols)
    !---------------------------------------------------------------------------
    ! Purpose: Set number of columns in the dynamics computational grid
    ! 
    ! Author: Jim Edwards
    !---------------------------------------------------------------------------
    use dimensions_mod, only: nelem, nelemd, np, npsq
    use element_mod,    only: index_t
    use dof_mod,        only: UniqueCoords, UniquePoints
    use physconst,      only: pi
    !------------------------------Arguments------------------------------------
    integer, intent(in) :: NumUniqueCols
    !----------------------------Local-Variables--------------------------------
    integer                :: ie, i, j, ii
    logical, save          :: local_coords_initialized = .false.
    real(r8)               :: areaw(np, np)
    type(index_t), pointer :: idx     
    !---------------------------------------------------------------------------
    if (.not. local_coords_initialized) then
      ngcols_d = NumUniqueCols
      ! Calculate the master mapping between element points and file order
      if (associated(fdofp_local)) then
        call endrun('set_horiz_grid_cnt_d: fdofp_local already defined')
      end if
      allocate(fdofp_local(npsq, nelemd))
      fdofp_local =0 
      do ie = 1, nelemd
        idx => elem(ie)%idxP
        do ii = 1, idx%NumUniquePts
          i = idx%ia(ii)
          j = idx%ja(ii)          
          fdofp_local((np*(j-1))+i,ie) = (idx%UniquePtoffset+ii-1)
        end do
      end do
      allocate(pelat_deg (np*np*nelemd))
      allocate(pelon_deg (np*np*nelemd))
      allocate(pearea    (np*np*nelemd))
      allocate(pemap     (np*np*nelemd))
      pemap = 0
      ! Now, fill in the appropriate values
      ii = 1
      do ie = 1, nelemd
        areaw = 1.0_r8 / elem(ie)%rspheremp(:,:)         
        pearea(ii:ii+npsq-1) = reshape(areaw, (/ np*np /))
        pemap(ii:ii+npsq-1) = fdofp_local(:,ie)
        do j = 1, np
          do i = 1, np
            pelat_deg(ii) = elem(ie)%spherep(i, j)%lat * rad2deg
            pelon_deg(ii) = elem(ie)%spherep(i, j)%lon * rad2deg
            ii = ii + 1
          end do
        end do
      end do
      local_coords_initialized = .true.
    else if (ngcols_d /= NumUniqueCols) then
      call endrun('set_horiz_grid_cnt_d: NumUniqueCols /= ngcols_d')
    end if

    return
  end subroutine set_horiz_grid_cnt_d
  !
  !=================================================================================================
  !
  subroutine get_horiz_grid_d(nxy,clat_d_out,clon_d_out,area_d_out, &
                              wght_d_out,lat_d_out,lon_d_out,cost_d_out)
    !--------------------------------------------------------------------------- 
    ! Purpose: Return latitude and longitude (in radians), column surface
    !          area (in radians squared) and surface integration weights
    !          for global column indices that will be passed to/from
    !          physics. Optionally also return estimated physics 
    !          computational cost per global column for use in load 
    !          balancing.
    ! 
    ! Author: Jim Edwards
    !------------------------------Arguments------------------------------------
    integer,  intent(in )                   :: nxy           ! array sizes
    real(r8), intent(out),         optional :: clat_d_out(:) ! column latitudes
    real(r8), intent(out),         optional :: clon_d_out(:) ! column longitudes
    real(r8), intent(out), target, optional :: area_d_out(:) ! column surface area
    real(r8), intent(out), target, optional :: wght_d_out(:) ! column integration weight
    real(r8), intent(out),         optional :: lat_d_out(:)  ! column degree latitudes
    real(r8), intent(out),         optional :: lon_d_out(:)  ! column degree longitudes
    real(r8), intent(out),        optional :: cost_d_out(:) ! column cost
    !----------------------------Local-Variables--------------------------------
    real(r8), pointer :: area_d(:)
    real(r8), allocatable :: lat_d_rad_temp(:)
    real(r8), allocatable :: lon_d_rad_temp(:)
    real(r8), allocatable :: lat_d_deg_temp(:)
    real(r8), allocatable :: lon_d_deg_temp(:)
    integer :: i
    !---------------------------------------------------------------------------
    ! Check that nxy is correct
    if (fv_nphys > 0) then
      if (nxy < fv_nphys*fv_nphys*nelem) then
         write(iulog, *) 'GET_HORIZ_GRID_D: arrays too small; Passed',     &
            nxy, ', needs to be at least', fv_nphys*fv_nphys*nelem
         call endrun('GET_HORIZ_GRID_D: arrays too small')
      end if
    else
      if (nxy < ngcols_d) then
        write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
             nxy,' < ',ngcols_d,' ) '
        call endrun('GET_HORIZ_GRID_D: arrays not large enough')
      end if
    end if ! fv_nphys > 0
    !---------------------------------------------------------------------------
    ! Area and/or integration weights
    if ( present(area_d_out) ) then
      if (size(area_d_out) /= nxy) then
        call endrun('bad area_d_out array size in dyn_grid')
      end if
      area_d => area_d_out
      call compute_global_area(area_d)
    else if ( present(wght_d_out) ) then
      if (size(wght_d_out) /= nxy) then
        call endrun('bad wght_d_out array size in dyn_grid')
      end if
      area_d => wght_d_out
      call compute_global_area(area_d)
    endif
    ! if both area_d_out or wght_d_out are present, then only wght_d_out was computed above
    if ( present(area_d_out) .and. present(wght_d_out) ) then
      wght_d_out(:) = area_d_out(:)
    end if
    !---------------------------------------------------------------------------
    ! Lat/Lon coordinates
    if (present(clon_d_out)) then
      if (size(clon_d_out).ne.nxy) call endrun('bad clon_d_out array size in dyn_grid')
    end if
    if ( present(clat_d_out) ) then
      if (size(clat_d_out).ne.nxy) call endrun('bad clat_d_out array size in dyn_grid')
    end if

    if ( present(clat_d_out) .or. present(clon_d_out) ) then
      
      allocate( lat_d_rad_temp(nxy) )
      allocate( lon_d_rad_temp(nxy) )
      allocate( lat_d_deg_temp(nxy) )
      allocate( lon_d_deg_temp(nxy) )

      call compute_global_coords(lat_d_rad_temp, lon_d_rad_temp, &
                                 lat_d_deg_temp, lon_d_deg_temp)

      if (present(clat_d_out)) clat_d_out = lat_d_rad_temp
      if (present(clon_d_out)) clon_d_out = lon_d_rad_temp
      if (present(lat_d_out))   lat_d_out = lat_d_deg_temp
      if (present(lon_d_out))   lon_d_out = lon_d_deg_temp

      deallocate( lat_d_rad_temp )
      deallocate( lon_d_rad_temp )
      deallocate( lat_d_deg_temp )
      deallocate( lon_d_deg_temp )

    end if
    !---------------------------------------------------------------------------
    ! just a placeholder until a mechanism for cost_d_out is implemented
    if (present(cost_d_out)) then
      if (size(cost_d_out) .ne. nxy) then
        call endrun('bad cost_d_out array size in dyn_grid')
      else
        cost_d_out(:) = 1.0_r8
      end if
    end if
    !---------------------------------------------------------------------------
    return
  end subroutine get_horiz_grid_d
  !
  !=================================================================================================
  !
  subroutine define_cam_grids()
    use cam_grid_support, only: horiz_coord_t, horiz_coord_create
    use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
    use gllfvremap_mod,   only: gfr_f_get_area, gfr_f_get_latlon
    !----------------------------Local-Variables--------------------------------
    character(len=16)            :: gridname, latname, lonname, ncolname, areaname
    type(horiz_coord_t), pointer :: lat_coord
    type(horiz_coord_t), pointer :: lon_coord
    integer(iMap),       pointer :: grid_map_d(:,:)
    integer(iMap),       pointer :: grid_map_p(:,:)
    integer                      :: ie, i, j, k, mapind ! Loop variables
    real(r8)                     :: area_scm(1), lat, lon
    integer                      :: ncols_p_lcl         ! local column count
    integer                      :: ncols_p_gbl         ! global column count
    integer(iMap),       pointer :: physgrid_map(:)
    real(r8),            pointer :: physgrid_area(:)
    real(r8),        allocatable :: physgrid_lat(:)
    real(r8),        allocatable :: physgrid_lon(:)
    !---------------------------------------------------------------------------
    ! Create grid for data on HOMME GLL nodes
    !---------------------------------------------------------------------------
    ! When the FV physs grid is used the GLL grid output 
    ! variables will use coordinates with the '_d' suffix.
    gridname = 'GLL'
    if (fv_nphys > 0) then
      latname  = 'lat_d'
      lonname  = 'lon_d'
      ncolname = 'ncol_d'
      areaname = 'area_d'
    else
      latname  = 'lat'
      lonname  = 'lon'
      ncolname = 'ncol'
      areaname = 'area'
    end if

    lat_coord => horiz_coord_create(trim(latname), trim(ncolname), ngcols_d,  &
                                    'latitude', 'degrees_north', 1,           &
                                    size(pelat_deg), pelat_deg, map=pemap)
    lon_coord => horiz_coord_create(trim(lonname), trim(ncolname), ngcols_d,  &
                                    'longitude', 'degrees_east', 1,           &
                                    size(pelon_deg), pelon_deg, map=pemap)

    ! Map for dynamics GLL grid
    allocate(grid_map_d(3, SIZE(fdofP_local)))
    grid_map_d = 0
    mapind = 1
    do j = LBOUND(fdofP_local, 2), UBOUND(fdofP_local, 2) ! should be 1, nelemd
      do i = LBOUND(fdofP_local, 1), UBOUND(fdofP_local, 1) ! should be 1,npsq
        grid_map_d(1, mapind) = i
        grid_map_d(2, mapind) = j
        grid_map_d(3, mapind) = fdofP_local(i, j)
        mapind = mapind + 1
      end do
    end do

    ! The native HOMME GLL grid
    call cam_grid_register(trim(gridname), dyn_decomp, lat_coord, lon_coord,  &
                           grid_map_d,block_indexed=.false., unstruct=.true.)
    if (.not.single_column .or. scm_multcols) then
      call cam_grid_attribute_register(trim(gridname), trim(areaname),   &
                                'gll grid areas', trim(ncolname), pearea, pemap)
    else
      ! if single column model, then this attribute has to be handled 
      ! by assigning just the SCM point. Else, the model will bomb out
      ! when writing the header information to history output
      area_scm(1) = 1.0_r8 / elem(1)%rspheremp(1,1)
      call cam_grid_attribute_register(trim(gridname), trim(areaname),   &
                                    'gll grid areas', trim(ncolname), area_scm)
    end if ! .not. single_column

    call cam_grid_attribute_register(trim(gridname), 'np', '', np)
    call cam_grid_attribute_register(trim(gridname), 'ne', '', ne)

    ! Coordinate values and maps are copied into the coordinate and attribute objects.
    ! Locally allocated storage is no longer needed.
    deallocate(pelat_deg)
    deallocate(pelon_deg)

    ! These can be nullified since the grid object has the reference
    nullify(grid_map_d) 
    nullify(pearea)
    nullify(pemap)

    !---------------------------------------------------------------------------
    ! Create grid object for physics grid on the dynamics decomposition
    !---------------------------------------------------------------------------
    ! Following CAM6-SE the 'physgrid_d' grid is created and used to load the 
    ! PHIS field (i.e. surface topo) on the FV physics grid and then 
    ! interpolated to the GLL grid. This ensures consistent treatment of SGH.
    if (fv_nphys>0) then

      gridname = 'physgrid_d'
      latname  = 'lat'
      lonname  = 'lon'
      ncolname = 'ncol'
      areaname = 'area'

      ncols_p_lcl = fv_nphys * fv_nphys * nelemd
      ncols_p_gbl = fv_nphys * fv_nphys * nelem
      allocate(physgrid_map(ncols_p_lcl))
      allocate(physgrid_lat(ncols_p_lcl))
      allocate(physgrid_lon(ncols_p_lcl))
      allocate(physgrid_area(ncols_p_lcl))

      ! copy local grid properties
      do ie = 1,nelemd
        k = 1
        do j = 1,fv_nphys
          do i = 1,fv_nphys
            mapind = k + (ie-1)*fv_nphys*fv_nphys
            physgrid_map(mapind) = k + (elem(ie)%GlobalId-1)*fv_nphys*fv_nphys
            physgrid_area(mapind)= gfr_f_get_area(ie, i, j)
            call gfr_f_get_latlon(ie, i, j, lat, lon)
            lat = lat*rad2deg
            lon = lon*rad2deg
            physgrid_lat(mapind) = lat
            physgrid_lon(mapind) = lon
            k = k + 1
          end do ! i
        end do ! j
      end do ! ie

      lat_coord => horiz_coord_create(trim(latname),trim(ncolname),ncols_p_gbl,&
                                      'latitude', 'degrees_north', 1,          &
                                      ncols_p_lcl,physgrid_lat,map=physgrid_map)

      lon_coord => horiz_coord_create(trim(lonname),trim(ncolname),ncols_p_gbl,&
                                      'longitude', 'degrees_east', 1,          &
                                      ncols_p_lcl,physgrid_lon,map=physgrid_map)

      ! Map for physics grid
      allocate(grid_map_p(3, ncols_p_lcl))
      grid_map_p = 0_iMap
      mapind = 1
      do j = 1, nelemd
        do i = 1, fv_nphys*fv_nphys
          grid_map_p(1,mapind) = i
          grid_map_p(2,mapind) = j
          grid_map_p(3,mapind) = physgrid_map(mapind)
          mapind = mapind + 1
        end do ! i
      end do ! j

      ! create physics grid object
      call cam_grid_register(trim(gridname), physgrid_d, lat_coord, lon_coord, &
                             grid_map_p, block_indexed=.false., unstruct=.true.)
      call cam_grid_attribute_register(trim(gridname), trim(areaname),         &
                                       'physics grid areas', trim(ncolname),   &
                                       physgrid_area, map=physgrid_map)
      call cam_grid_attribute_register(trim(gridname), 'fv_nphys', '', fv_nphys)
      call cam_grid_attribute_register(trim(gridname), 'ne',       '', ne)

      deallocate(physgrid_lat)
      deallocate(physgrid_lon)
      nullify(physgrid_area)
      nullify(physgrid_map)
      nullify(grid_map_p)

    end if ! fv_nphys>0
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    nullify(lat_coord) ! Belongs to grid
    nullify(lon_coord) ! Belongs to grid
  end subroutine define_cam_grids
  !
  !=================================================================================================
  !
  subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)
    use cam_grid_support, only: max_hcoordname_len
    !------------------------------Arguments------------------------------------
    character(len=max_hcoordname_len),          intent(out) :: gridname
    character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)
    !---------------------------------------------------------------------------
    if (fv_nphys > 0) then
      gridname = 'physgrid_d'
      allocate(grid_attribute_names(3))
      grid_attribute_names(1) = 'area'
      grid_attribute_names(2) = 'fv_nphys'
      grid_attribute_names(3) = 'ne'
    else
      gridname = 'GLL'
      allocate(grid_attribute_names(3))
      grid_attribute_names(1) = 'area'
      grid_attribute_names(2) = 'np'
      grid_attribute_names(3) = 'ne'
    end if ! fv_nphys > 0

  end subroutine physgrid_copy_attributes_d
  !
  !=================================================================================================
  !
  subroutine gblocks_init()
    !----------------------------Local-Variables--------------------------------
    integer :: ie, p
    integer :: ibuf
    integer :: ierr
    integer :: rdispls(npes)
    integer :: recvcounts(npes)
    integer :: gid(npes)
    integer :: lid(npes)
    !---------------------------------------------------------------------------
    if (.not.allocated(gblocks)) then
      if (masterproc) then
        write(iulog, *) 'INFO: Non-scalable action: Allocating global blocks in SE dycore.'
      end if
      allocate(gblocks(nelem))
      do ie=1,nelem
        gblocks(ie)%Owner          = -1
        gblocks(ie)%UniquePtOffset = -1
        gblocks(ie)%NumUniqueP     = -1
        gblocks(ie)%LocalID        = -1
      end do
    end if

    do ie = 1,nelemdmax        
      if (ie<=nelemd) then
        rdispls(iam+1)=elem(ie)%idxP%UniquePtOffset-1
        gid(iam+1)=elem(ie)%GlobalID
        lid(iam+1)=ie
        recvcounts(iam+1)=elem(ie)%idxP%NumUniquePts
      else     
        rdispls(iam+1) = 0
        recvcounts(iam+1) = 0
        gid(iam+1)=0
      end if

      ibuf=lid(iam+1)
      call mpi_allgather(ibuf, 1, mpi_integer, lid       , 1, mpi_integer, mpicom, ierr)
      ibuf=gid(iam+1)
      call mpi_allgather(ibuf, 1, mpi_integer, gid       , 1, mpi_integer, mpicom, ierr)
      ibuf=rdispls(iam+1)
      call mpi_allgather(ibuf, 1, mpi_integer, rdispls   , 1, mpi_integer, mpicom, ierr)
      ibuf=recvcounts(iam+1)
      call mpi_allgather(ibuf, 1, mpi_integer, recvcounts, 1, mpi_integer, mpicom, ierr)

      do p=1,npes
        if(gid(p)>0) then
          gblocks(gid(p))%UniquePtOffset=rdispls(p)+1
          gblocks(gid(p))%NumUniqueP=recvcounts(p)
          gblocks(gid(p))%LocalID=lid(p)
          gblocks(gid(p))%Owner=p-1
        end if
      end do ! p
    end do ! ie
    gblocks_need_initialized=.false.

  end subroutine gblocks_init
  !
  !=================================================================================================
  !
  subroutine compute_global_area(area_d)
    use dof_mod,                only: UniqueCoords, UniquePoints
    use gllfvremap_mod,         only: gfr_f_get_area
    !------------------------------Arguments------------------------------------
    real(r8), pointer :: area_d(:)
    !----------------------------Local-Variables--------------------------------
    real(r8), dimension(np,np)  :: areaw
    real(r8), allocatable       :: area_fv(:,:)
    real(r8), allocatable       :: rbuf(:)
    integer,  dimension(npes)   :: displace  ! MPI data displacement for gathering
    integer,  dimension(npes)   :: recvcnts  ! MPI send buffer count for gathering
    integer,  allocatable       :: col_id_global(:)
    integer,  allocatable       :: col_id_local(:)
    real(r8), allocatable       :: area_local(:)
    integer  :: ncol_fv_gbl, ncol_fv_lcl
    integer  :: ie, sb, eb, i, j, ip, fv_cnt, icol
    integer  :: ierr
    integer  :: ibuf
    !---------------------------------------------------------------------------
    if (masterproc) then
      write(iulog,*) 'INFO: Non-scalable action: Computing global area in SE dycore.'
    end if

    if (fv_nphys > 0) then
      
      ncol_fv_gbl = fv_nphys*fv_nphys*nelem
      ncol_fv_lcl = fv_nphys*fv_nphys*nelemd
      allocate(rbuf(ncol_fv_gbl))
      allocate(col_id_local(ncol_fv_lcl))
      allocate(col_id_global(ncol_fv_gbl))
      allocate(area_local(ncol_fv_gbl))

      ! Get area for local blocks
      icol = 1
      do ie = 1,nelemd
        do j = 1,fv_nphys
          do i = 1,fv_nphys
            area_local(icol) = gfr_f_get_area(ie, i, j)
            icol = icol+1
          end do ! i
        end do ! j
      end do ! ie

      ! gather send buffer count as local cell count
      recvcnts(:) = 0
      call mpi_allgather(ncol_fv_lcl, 1, mpi_integer, recvcnts, 1, mpi_integer, mpicom, ierr)

      ! determine displacement for MPI gather
      displace(1) = 0
      do ip = 2,npes
        displace(ip) = displace(ip-1) + recvcnts(ip-1)
      end do
      ! Check to make sure we counted correctly
      if (masterproc) then
        if ( displace(npes) + recvcnts(npes) /= ncol_fv_gbl ) then
          call endrun('compute_global_area: bad MPI displace array size')
        end if
      end if ! masterproc

      ! gather element IDs for sorting
      col_id_global(:) = -1
      icol = 1
      do ie = 1,nelemd
        fv_cnt = 1
        do j = 1,fv_nphys
          do i = 1,fv_nphys
            col_id_local(icol) = (elem(ie)%GlobalId-1)*fv_nphys*fv_nphys + fv_cnt
            fv_cnt = fv_cnt+1
            icol = icol+1
          end do
        end do
      end do
      call mpi_allgatherv( col_id_local(1:ncol_fv_lcl), recvcnts(iam+1), mpi_integer, col_id_global, &
                           recvcnts(:), displace(:), mpi_integer, mpicom, ierr)


      call mpi_allgatherv( area_local(:), recvcnts(iam+1), mpi_real8, rbuf, &
                           recvcnts(:), displace(:), mpi_real8, mpicom, ierr)

      ! sort according to global element ID
      do icol = 1,ncol_fv_gbl
        area_d( col_id_global(icol) ) = rbuf(icol)
      end do

      deallocate(rbuf)
      deallocate(area_local)
      deallocate(col_id_local)
      deallocate(col_id_global)

    else ! physics is on GLL grid
    
      allocate(rbuf(ngcols_d))
      do ie = 1, nelemdmax        
        if(ie <= nelemd) then
          displace(iam+1) = elem(ie)%idxp%UniquePtOffset-1
          recvcnts(iam+1) = elem(ie)%idxP%NumUniquePts
          eb = displace(iam+1) + elem(ie)%idxp%NumUniquePts
          areaw = 1.0_r8 / elem(ie)%rspheremp(:,:)         
          call UniquePoints(elem(ie)%idxP, areaw, area_d(displace(iam+1)+1:eb))
        else
          displace(iam+1) = 0
          recvcnts(iam+1) = 0
        end if
        ibuf = displace(iam+1)
        call mpi_allgather(ibuf, 1, mpi_integer, displace, 1, mpi_integer, mpicom, ierr)
        ibuf = recvcnts(iam+1)
        call mpi_allgather(ibuf, 1, mpi_integer, recvcnts, 1, mpi_integer, mpicom, ierr)
        sb = displace(iam+1) + 1
        eb = displace(iam+1) + recvcnts(iam+1)
        rbuf(1:recvcnts(iam+1)) = area_d(sb:eb)
        call mpi_allgatherv(rbuf, recvcnts(iam+1), mpi_real8, area_d,       &
                            recvcnts(:), displace(:), mpi_real8, mpicom, ierr)
      end do ! ie
      deallocate(rbuf)

    end if ! fv_nphys > 0

  end subroutine compute_global_area
  !
  !=================================================================================================
  !
  subroutine compute_global_coords(clat, clon, lat_out, lon_out, corner_lat_out, corner_lon_out)
    use dof_mod,                only: UniqueCoords, UniquePoints
    use gllfvremap_mod,         only: gfr_f_get_latlon, gfr_f_get_corner_latlon
    !------------------------------Arguments------------------------------------
    real(r8),           intent(out) :: clat(:)                ! radians
    real(r8),           intent(out) :: clon(:)                ! radians
    real(r8), optional, intent(out) :: lat_out(:)             ! degrees
    real(r8), optional, intent(out) :: lon_out(:)             ! degrees
    real(r8), optional, intent(out) :: corner_lat_out(:,:)    ! degrees
    real(r8), optional, intent(out) :: corner_lon_out(:,:)    ! degrees
    !----------------------------Local-Variables--------------------------------
    real(r8), allocatable     :: rbuf(:)
    integer,  dimension(npes) :: displace  ! MPI data displacement for gathering
    integer,  dimension(npes) :: recvcnts  ! MPI send buffer count for gathering
    integer,  allocatable     :: col_id_global(:)
    integer,  allocatable     :: col_id_local(:)
    real(r8), allocatable     :: lat_rad_local(:)
    real(r8), allocatable     :: lon_rad_local(:)
    real(r8), allocatable     :: corner_lat_rad_local(:,:)
    real(r8), allocatable     :: corner_lon_rad_local(:,:)
    integer  :: ncol_fv_gbl, ncol_fv_lcl
    integer  :: ie, sb, eb, j, i, ip, fv_cnt, icol, c
    integer  :: ierr
    integer  :: ibuf
    !---------------------------------------------------------------------------
    if (masterproc) then
      write(iulog,*) 'INFO: Non-scalable action: Computing global coords in SE dycore.'
    end if

    if (fv_nphys > 0) then

      ncol_fv_gbl = fv_nphys*fv_nphys*nelem
      ncol_fv_lcl = fv_nphys*fv_nphys*nelemd
      allocate(rbuf(ncol_fv_gbl))
      allocate(col_id_local(ncol_fv_lcl))
      allocate(col_id_global(ncol_fv_gbl))
      allocate(lat_rad_local(ncol_fv_lcl))
      allocate(lon_rad_local(ncol_fv_lcl))

      ! calculate coordinates for local blocks
      icol = 1
      do ie = 1,nelemd
        do j = 1,fv_nphys
          do i = 1,fv_nphys
            call gfr_f_get_latlon(ie, i, j, lat_rad_local(icol), lon_rad_local(icol))
            icol = icol+1
          end do ! i
        end do ! j
      end do ! ie

      ! gather send buffer count as local cell count
      recvcnts(:) = 0
      call mpi_allgather(ncol_fv_lcl, 1, mpi_integer, recvcnts, 1, mpi_integer, mpicom, ierr)

      ! determine displacement for MPI gather
      displace(1) = 0
      do ip = 2,npes
        displace(ip) = displace(ip-1) + recvcnts(ip-1)
      end do

      ! Check to make sure we counted correctly
      if (masterproc) then
        if ( displace(npes) + recvcnts(npes) /= ncol_fv_gbl ) then
          call endrun('compute_global_coords: bad MPI displace array size')
        end if
      end if

      ! gather element IDs for sorting
      col_id_global(:) = -1
      icol = 1
      do ie = 1,nelemd
        fv_cnt = 1
        do j = 1,fv_nphys
          do i = 1,fv_nphys
            col_id_local(icol) = (elem(ie)%GlobalId-1)*fv_nphys*fv_nphys + fv_cnt
            fv_cnt = fv_cnt+1
            icol = icol+1
          end do
        end do
      end do
      call mpi_allgatherv( col_id_local(1:ncol_fv_lcl), recvcnts(iam+1), mpi_integer, col_id_global, &
                           recvcnts(:), displace(:), mpi_integer, mpicom, ierr)

      ! Gather global latitudes
      call mpi_allgatherv( lat_rad_local(:), recvcnts(iam+1), mpi_real8, rbuf, &
                           recvcnts(:), displace(:), mpi_real8, mpicom, ierr)
      
      ! sort latitude according to global element ID
      do icol = 1,ncol_fv_gbl
        clat( col_id_global(icol) ) = rbuf(icol)
      end do

      ! Gather global longitudes
      call mpi_allgatherv( lon_rad_local(:), recvcnts(iam+1), mpi_real8, rbuf, &
                           recvcnts(:), displace(:), mpi_real8, mpicom, ierr)
      
      ! sort longitude according to global element ID
      do icol = 1,ncol_fv_gbl
        clon( col_id_global(icol) ) = rbuf(icol)
      end do

      ! Create version in degrees if requested
      if (present(lat_out)) lat_out(:) = clat(:) * rad2deg
      if (present(lon_out)) lon_out(:) = clon(:) * rad2deg

      !----------------------------------------------------
      ! Get cell corners (only needed for writing scrip file)
      !----------------------------------------------------
      if ( present(corner_lat_out) .or. present(corner_lon_out) ) then
        allocate(corner_lat_rad_local(ncol_fv_lcl,4))
        allocate(corner_lon_rad_local(ncol_fv_lcl,4))
        icol = 1
        do ie = 1,nelemd
          do j = 1,fv_nphys
            do i = 1,fv_nphys
              do c = 1,4
                call gfr_f_get_corner_latlon(ie, i, j, c, &
                     corner_lat_rad_local(icol,c), corner_lon_rad_local(icol,c))
              end do
              icol = icol+1
            end do ! i
          end do ! j
        end do ! ie
        ! Gather corner coordinates (one corner at a time)
        do c = 1,4
          ! Gather latitude
          call mpi_allgatherv( corner_lat_rad_local(:,c), recvcnts(iam+1), mpi_real8, rbuf, &
                               recvcnts(:), displace(:), mpi_real8, mpicom, ierr)
          ! sort latitude 
          do icol = 1,ncol_fv_gbl
            corner_lat_out( col_id_global(icol), c) = rbuf(icol) * rad2deg
          end do
          ! Gather longitude
          call mpi_allgatherv( corner_lon_rad_local(:,c), recvcnts(iam+1), mpi_real8, rbuf, &
                               recvcnts(:), displace(:), mpi_real8, mpicom, ierr)
          ! sort longitude 
          do icol = 1,ncol_fv_gbl
            corner_lon_out( col_id_global(icol), c) = rbuf(icol) * rad2deg
          end do
        end do ! c
        ! Deallocate stuff for corners
        deallocate(corner_lat_rad_local)
        deallocate(corner_lon_rad_local)
      end if ! present(corner_lat_out) .or. present(corner_lon_out)
      !----------------------------------------------------
      !----------------------------------------------------

      ! Deallocate stuff
      deallocate(rbuf)
      deallocate(lat_rad_local)
      deallocate(lon_rad_local)
      deallocate(col_id_local)
      deallocate(col_id_global)

    else ! physics is on GLL grid

      allocate(rbuf(ngcols_d))

      clat(:) = -iam
      clon(:) = -iam
      if (present(lon_out)) then
        lon_out(:) = -iam
      end if
      if (present(lat_out)) then
        lat_out(:) = -iam
      end if

      do ie = 1, nelemdmax
        if(ie <= nelemd) then
          displace(iam+1) = elem(ie)%idxp%UniquePtOffset-1
          eb = displace(iam+1) + elem(ie)%idxp%NumUniquePts
          recvcnts(iam+1) = elem(ie)%idxP%NumUniquePts
          call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep,                    &
               clat(displace(iam+1)+1:eb),                                       &
               clon(displace(iam+1)+1:eb))
          if (present(lat_out)) then
            lat_out(displace(iam+1)+1:eb) = clat(displace(iam+1)+1:eb) * rad2deg
          end if
          if (present(lon_out)) then
            lon_out(displace(iam+1)+1:eb) = clon(displace(iam+1)+1:eb) * rad2deg
          end if
        else
          displace(iam+1) = 0
          recvcnts(iam+1) = 0
        end if
        ibuf = displace(iam+1)
        call mpi_allgather(ibuf, 1, mpi_integer, displace, &
             1, mpi_integer, mpicom, ierr)

        ibuf = recvcnts(iam+1)
        call mpi_allgather(ibuf, 1, mpi_integer, recvcnts, &
             1, mpi_integer, mpicom, ierr)

        sb = displace(iam+1) + 1
        eb = displace(iam+1) + recvcnts(iam+1)

        rbuf(1:recvcnts(iam+1)) = clat(sb:eb)
        call mpi_allgatherv(rbuf, recvcnts(iam+1), mpi_real8, clat,            &
               recvcnts(:), displace(:), mpi_real8, mpicom, ierr)
        if (present(lat_out)) then
          rbuf(1:recvcnts(iam+1)) = lat_out(sb:eb) 
          call mpi_allgatherv(rbuf, recvcnts(iam+1), mpi_real8, lat_out,       &
               recvcnts(:), displace(:), mpi_real8, mpicom, ierr)
        end if

        rbuf(1:recvcnts(iam+1)) = clon(sb:eb)
        call mpi_allgatherv(rbuf, recvcnts(iam+1), mpi_real8, clon,            &
               recvcnts(:), displace(:), mpi_real8, mpicom, ierr)
        if (present(lon_out)) then
          rbuf(1:recvcnts(iam+1)) = lon_out(sb:eb) 
          call mpi_allgatherv(rbuf, recvcnts(iam+1), mpi_real8, lon_out,       &
               recvcnts(:), displace(:), mpi_real8, mpicom, ierr)
        end if
      end do
      
      deallocate(rbuf)

    end if ! fv_nphys > 0

  end subroutine compute_global_coords
  !
  !=================================================================================================
  !
  function get_dyn_grid_parm_real2d(name) result(rval)
    !------------------------------Arguments------------------------------------
    character(len=*), intent(in) :: name
    real(r8),            pointer :: rval(:,:)
    !---------------------------------------------------------------------------
    if(name.eq.'clon') then
      call endrun('get_dyn_grid_parm_real2d: clon not supported, use get_horiz_grid_d')
    else if(name.eq.'londeg') then
      ! This is never set so I'm calling bogus and stomping calls --goldy
      call endrun('get_dyn_grid_parm_real2d: londeg not defined')
    else
      nullify(rval)
    end if
  end function get_dyn_grid_parm_real2d
  !
  !=================================================================================================
  !
  function get_dyn_grid_parm_real1d(name) result(rval)
    !------------------------------Arguments------------------------------------
    character(len=*), intent(in) :: name
    real(r8),            pointer :: rval(:)
    !---------------------------------------------------------------------------
    if(name.eq.'w') then
      if (.not. associated(w)) then
        call endrun('get_dyn_grid_parm_real1d: w not defined')
      end if
      rval => w
    else if(name.eq.'clat') then
      call endrun('get_dyn_grid_parm_real1d: clat not supported, use get_horiz_grid_d')
    else if(name.eq.'latdeg') then
      ! This is never set so I'm calling bogus and stomping calls --goldy
      call endrun('get_dyn_grid_parm_real1d: latdeg not defined')
    else
      nullify(rval)
    end if
  end function get_dyn_grid_parm_real1d
  !
  !=================================================================================================
  !
  integer function get_dyn_grid_parm(name) result(ival)
    use pmgrid,          only: beglat, endlat, plat, plon, plev, plevp
    use interpolate_mod, only: get_interp_parameter
    !------------------------------Arguments------------------------------------
    character(len=*), intent(in) :: name
    !---------------------------------------------------------------------------
    if(name.eq.'ne') then
      ival = ne
    else if(name.eq.'np') then
      ival = np
    else if(name.eq.'npsq') then
      ival = npsq
    else if(name.eq.'nelemd') then
      ival = nelemd
    else if(name.eq.'beglat') then
      ival = beglat
    else if(name.eq.'endlat') then
      ival = endlat
      ! The following four are required for consistancy in cam_history
    else if(name.eq.'beglonxy') then
      ival = 1
    else if(name.eq.'endlonxy') then
      ival = npsq
    else if(name.eq.'beglatxy') then
      ival=1
    else if(name.eq.'endlatxy') then
      ival=nelemd
    else if(name.eq.'plat') then
      ival = plat
    else if(name.eq.'plon') then
      if (fv_nphys>0) then
        ival = fv_nphys*fv_nphys*nelem
      else
        ival = ngcols_d
      end if 
    else if(name.eq.'plev') then
      ival = plev
    else if(name.eq.'plevp') then
      ival = plevp
    else if(name.eq.'nlon') then
      ival = get_interp_parameter('nlon')
    else if(name.eq.'nlat') then
      ival = get_interp_parameter('nlat')
    else if(name.eq.'num_grids') then
      ival = 2
    else
      ival = -1
    end if
  end function get_dyn_grid_parm
  !
  !=================================================================================================
  !
  subroutine dyn_grid_get_pref(pref_edge, pref_mid, num_pr_lev)
    !---------------------------------------------------------------------------
    ! Purpose: return reference pressures for the dynamics grid
    !---------------------------------------------------------------------------
    use pmgrid, only: plev
    use hycoef, only: hypi, hypm, nprlev
    !------------------------------Arguments--------------------------------
    real(r8), intent(out) :: pref_edge(:) ! reference pressure at layer edges (Pa)
    real(r8), intent(out) :: pref_mid(:)  ! reference pressure at layer midpoints (Pa)
    integer,  intent(out) :: num_pr_lev   ! number of top levels using pure pressure representation
    !----------------------------Local-Variables--------------------------------
    integer :: k
    !---------------------------------------------------------------------------
    do k = 1, plev
      pref_edge(k) = hypi(k)
      pref_mid(k)  = hypm(k)
    end do
    pref_edge(plev+1) = hypi(plev+1)

    num_pr_lev = nprlev

  end subroutine dyn_grid_get_pref
  !
  !=================================================================================================
  !
  subroutine dyn_grid_find_gcols( lat, lon, nclosest, owners, col, lbk, rlat, rlon, idyn_dists ) 
    !---------------------------------------------------------------------------
    ! This returns the lat/lon information (and corresponding MPI task 
    ! numbers (owners))  of the global model grid columns nearest to the 
    ! input satellite coordinate (lat,lon)
    !---------------------------------------------------------------------------
    use shr_const_mod,  only: SHR_CONST_REARTH
    !------------------------------Arguments------------------------------------
    real(r8),          intent(in ) :: lat
    real(r8),          intent(in ) :: lon
    integer,           intent(in ) :: nclosest
    integer,           intent(out) :: owners(nclosest)
    integer,           intent(out) :: col(nclosest)
    integer,           intent(out) :: lbk(nclosest)
    real(r8),optional, intent(out) :: rlon(nclosest)
    real(r8),optional, intent(out) :: rlat(nclosest)
    real(r8),optional, intent(out) :: idyn_dists(nclosest)
    !----------------------------Local-Variables--------------------------------
    real(r8) :: dist            ! distance (in radians**2 from lat, lon)
    real(r8) :: latr            ! lat inputs converted to radians
    real(r8) :: lonr            ! lon inputs converted to radians
    integer  :: blockid(1)
    integer  :: bcid(1)
    integer  :: lclblockid(1)
    integer  :: i,j,k, ii
    real(r8), allocatable :: clat_d(:)
    real(r8), allocatable :: clon_d(:)
    real(r8), allocatable :: distmin(:)
    integer,  allocatable :: igcol(:)
    !---------------------------------------------------------------------------
    latr = lat/rad2deg
    lonr = lon/rad2deg

    allocate( clat_d(1:ngcols_d) )
    allocate( clon_d(1:ngcols_d) )
    allocate( igcol(nclosest) )
    allocate( distmin(nclosest) )

    call get_horiz_grid_d(ngcols_d, clat_d_out=clat_d, clon_d_out=clon_d)

    igcol(:)    = -999
    distmin(:) = 1.e10_r8

    do i = 1,ngcols_d

      ! Use the Spherical Law of Cosines to find the great-circle distance.
      dist = acos(sin(latr) * sin(clat_d(i)) + cos(latr) * cos(clat_d(i)) * cos(clon_d(i) - lonr)) * SHR_CONST_REARTH
      do j = nclosest, 1, -1
        if (dist < distmin(j)) then

          if (j < nclosest) then
            distmin(j+1) = distmin(j)
            igcol(j+1)    = igcol(j)
          end if

          distmin(j) = dist
          igcol(j)    = i
        else
          exit
        end if
      enddo

    enddo

    do i = 1,nclosest

      call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
      owners(i) = get_block_owner_d(blockid(1))

      if (owners(i)==iam) then
        lbk(i) = lclblockid(1)
        ii = igcol(i)-elem(lbk(i))%idxp%UniquePtoffset+1
        k=elem(lbk(i))%idxp%ia(ii)
        j=elem(lbk(i))%idxp%ja(ii)
        col(i) = k+(j-1)*np
      else
        lbk(i) = -1
        col(i) = -1
      endif
      if ( present(rlat) ) rlat(i) = clat_d(igcol(i)) * rad2deg
      if ( present(rlon) ) rlon(i) = clon_d(igcol(i)) * rad2deg

      if (present(idyn_dists)) then
        idyn_dists(i) = distmin(i)
      end if

    end do

    deallocate( clat_d )
    deallocate( clon_d )
    deallocate( igcol )
    deallocate( distmin )

  end subroutine dyn_grid_find_gcols
  !
  !=================================================================================================
  !
  subroutine dyn_grid_get_colndx( igcol, nclosest, owners, col, lbk ) 
    !------------------------------Arguments------------------------------------
    integer, intent(in)  :: nclosest
    integer, intent(in)  :: igcol(nclosest)
    integer, intent(out) :: owners(nclosest)
    integer, intent(out) :: col(nclosest)
    integer, intent(out) :: lbk(nclosest)
    !----------------------------Local-Variables--------------------------------
    integer  :: i,j,k, ii
    integer  :: blockid(1)
    integer  :: bcid(1)
    integer  :: lclblockid(1)
    !---------------------------------------------------------------------------

    if (fv_nphys > 0) then
      call endrun('dyn_grid_get_colndx: not implemented for the FV physics grid')
    end if

    do i = 1,nclosest

       call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
       owners(i) = get_block_owner_d(blockid(1))
       
       if (owners(i)==iam) then
          lbk(i) = lclblockid(1)
          ii = igcol(i)-elem(lbk(i))%idxp%UniquePtoffset+1
          k=elem(lbk(i))%idxp%ia(ii)
          j=elem(lbk(i))%idxp%ja(ii)
          col(i) = k+(j-1)*np
       else
          lbk(i) = -1
          col(i) = -1
       endif
       
    end do

  end subroutine dyn_grid_get_colndx
  !
  !=================================================================================================
  !
  subroutine dyn_grid_get_elem_coords( ie, rlon, rlat, cdex )
    !---------------------------------------------------------------------------
    ! Purpose: this returns coordinates of a specified block element of the dyn grid
    !---------------------------------------------------------------------------
    use dof_mod, only: UniqueCoords
    !------------------------------Arguments------------------------------------
    integer,           intent(in ) :: ie      ! block element index
    real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the element
    real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the element
    integer, optional, intent(out) :: cdex(:) ! global column index
    !----------------------------Local-Variables--------------------------------
    integer :: sb,eb, ii, i,j, icol, igcol
    real(r8), allocatable :: clat(:)
    real(r8), allocatable :: clon(:)
    !---------------------------------------------------------------------------
    if (fv_nphys > 0) then
      call endrun('dyn_grid_get_elem_coords: not implemented for the FV physics grid')
    end if

    sb = elem(ie)%idxp%UniquePtOffset
    eb = sb + elem(ie)%idxp%NumUniquePts-1

    allocate( clat(sb:eb), clon(sb:eb) )
    call UniqueCoords( elem(ie)%idxP, elem(ie)%spherep, clat(sb:eb), clon(sb:eb) )

    if (present(cdex)) cdex(:) = -1
    if (present(rlat)) rlat(:) = -999._r8
    if (present(rlon)) rlon(:) = -999._r8

    do ii=1,elem(ie)%idxp%NumUniquePts
      i=elem(ie)%idxp%ia(ii)
      j=elem(ie)%idxp%ja(ii)
      icol = i+(j-1)*np
      igcol = elem(ie)%idxp%UniquePtoffset+ii-1
      if (present(cdex)) cdex(icol) = igcol
      if (present(rlat)) rlat(icol) = clat( igcol )
      if (present(rlon)) rlon(icol) = clon( igcol )
    end do

    deallocate( clat, clon )

  end subroutine dyn_grid_get_elem_coords
  !
  !=================================================================================================
  !
  subroutine fv_physgrid_init()
    use gllfvremap_mod,         only: gfr_init

    call gfr_init(par, elem, fv_nphys)
  end subroutine fv_physgrid_init
  !
  !=================================================================================================
  !
  subroutine fv_physgrid_final()
    use gllfvremap_mod,   only: gfr_finish

    call gfr_finish()
  end subroutine fv_physgrid_final
  !
  !=================================================================================================
  !
end module dyn_grid
