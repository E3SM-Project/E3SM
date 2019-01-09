module q3d_coupling

!-------------------------------------------------------------------------------
! Q3D vGCM - physics coupling module
!-------------------------------------------------------------------------------

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use ppgrid,         only: pver
  use constituents,   only: pcnst

  use dyn_grid,       only: TimeLevel
  use dp_mapping,     only: nphys_pts

  use cam_logfile,    only: iulog
  use shr_sys_mod,    only: shr_sys_flush
  use perf_mod,       only: t_startf, t_stopf, t_barrierf
  use cam_abortutils, only: endrun

  use parallel_mod,   only: par
  use hybrid_mod,     only: config_thread_region, get_loop_ranges, hybrid_t
  use dimensions_mod, only: nelemd, nlev, nc, qsize, ntrac, fv_nphys

  use parmsld,        only: nhalo_vGCM

  implicit none
  private
  save

  public :: q3d_set_indexing ! One-time index setup
  public :: q3d_halo_data    ! Called to retrieve state information
  public :: q3d_halo_coord   ! Called to retrieve coordinate information

  ! Indices for field arrays returned from q3d_halo_data
  ! Midpoint fields (nlev)
  integer, parameter, public :: halo_index_u_con   = 1
  integer, parameter, public :: halo_index_v_con   = 2
  integer, parameter, public :: halo_index_Temp    = 3
  integer, parameter, public :: halo_index_OMEGA   = 4
  integer, parameter, public :: halo_index_DPDRY   = 5
  integer, parameter, public :: halo_index_Q       = 6
  integer, parameter, public :: halo_num_midpoint  = halo_index_Q + pcnst -1
  ! Interface fields (nlev+1)
  integer, parameter, public :: halo_index_pint    = 1
  integer, parameter, public :: halo_index_zi      = 2
  integer, parameter, public :: halo_num_interface = halo_index_zi
  ! Surface fields
  integer, parameter, public :: halo_index_PHIS    = 1
  integer, parameter, public :: halo_num_surface   = halo_index_PHIS
  integer, parameter, public :: halo_index_lon     = 1
  integer, parameter, public :: halo_index_lat     = 2

  type, public :: halo_index_t
     integer :: gcid       ! Source column global index
     integer :: ghost_id   ! Ghost column global index
     integer :: ie         ! Source column local element index
     integer :: ioff       ! Source column local i offset
     integer :: joff       ! Source column local j offset
     integer :: ighost     ! i offset from source column for ghost column
     integer :: jghost     ! j offset from source column for ghost column
     integer :: send_c1_indices(-nhalo_vGCM:nhalo_vGCM) ! Indices of vGCM sending column(s)
     integer :: send_c2_indices(-nhalo_vGCM:nhalo_vGCM) ! Indices of vGCM sending column(s)
  end type halo_index_t

  type, public :: element_index_t
     integer :: gcid       ! Source column global index
     integer :: ie         ! Source column local element index
     integer :: ioff       ! Source column local i offset
     integer :: joff       ! Source column local j offset
  end type element_index_t

  interface q3d_set_indexing
     module procedure q3d_set_indexing_ghost
     module procedure q3d_set_indexing_element
  end interface q3d_set_indexing

!============================================================================
CONTAINS
!============================================================================

   subroutine q3d_halo_data(elem_fill_index, ghost_fill_index, phys_outs, phys_outm, phys_outi)

     use spmd_dyn,               only: local_dp_map
     use spmd_utils,             only: iam
     use hycoef,                 only: hyai, ps0
     use fvm_control_volume_mod, only: n0_fvm
     use fvm_mapping,            only: dyn2phys_vector, dyn2phys_all_vars
     use time_mod,               only: timelevel_qdp
     use control_mod,            only: qsplit
     use dimensions_mod,         only: nhc_phys,nhr_phys
     use dyn_grid,               only: elem, fvm

     ! arguments
     type(element_index_t), intent(in)    :: elem_fill_index(:)
     type(halo_index_t),    intent(in)    :: ghost_fill_index(:)
     real(kind=r8),         intent(inout) :: phys_outs(:,:)   ! Surface
     real(kind=r8),         intent(inout) :: phys_outm(:,:,:) ! Midpoint
     real(kind=r8),         intent(inout) :: phys_outi(:,:,:) ! Interface

     ! LOCAL VARIABLES
     integer                      :: ie               ! indices over elements
     real (kind=r8),  allocatable :: ps_tmp(:,:)      ! temp array to hold ps
     real (kind=r8),  allocatable :: dp3d_tmp(:,:,:)  ! temp array to hold dp3d
     real (kind=r8),  allocatable :: phis_tmp(:,:)    ! temp array to hold phis
     real (kind=r8),  allocatable :: T_tmp(:,:,:)     ! temp array to hold T
     real (kind=r8),  allocatable :: uv_tmp(:,:,:,:)  ! temp array to hold wind
     real (kind=r8),  allocatable :: q_tmp(:,:,:,:)   ! temp for advected constituents
     real (kind=r8),  allocatable :: omega_tmp(:,:,:) ! temp array to hold omega
     real (kind=r8),  allocatable :: fld_phys(:,:,:,:,:),fld_phys_int(:,:,:,:,:),fld_phys_2d(:,:,:,:,:)
     integer                      :: i, j, k, m, index, ioff

     integer                      :: tl_f, tl_qdp_np0, tl_qdp_np1
     integer                      :: nets, nete
     type(hybrid_t)               :: hybrid

     logical                      :: positive_definite(halo_num_midpoint)

     character(len=*),  parameter :: sub = 'q3d_halo_data'
!----------------------------------------------------------------------------

     positive_definite = .false.
     positive_definite(halo_index_Q:) = .true.

     tl_f = TimeLevel%n0
     call TimeLevel_Qdp(TimeLevel, qsplit, tl_qdp_np0,tl_qdp_np1)

     if (fv_nphys /= 3) then
        call endrun(sub//': Q3D requires fv_nphys = 3')
     end if
     if (.not. local_dp_map) then
        call endrun(sub//': Q3D only supports local maps')
     end if

     ! Allocate temporary arrays to hold data for physics decomposition
     allocate(ps_tmp(nphys_pts,nelemd))
     allocate(dp3d_tmp(nphys_pts,pver,nelemd))
     allocate(phis_tmp(nphys_pts,nelemd))
     allocate(T_tmp(nphys_pts,pver,nelemd))
     allocate(uv_tmp(nphys_pts,2,pver,nelemd))
     allocate(q_tmp(nphys_pts,pver,pcnst,nelemd))
     allocate(omega_tmp(nphys_pts,pver,nelemd))

     call t_startf('dyn2phys')
     if (iam < par%nprocs) then

        !******************************************************************
        ! physics runs on an FVM grid: map GLL vars to physics grid
       !******************************************************************      
       do ie = 1, nelemd
           ! note that the fvm halo has been filled in prim_run_subcycle
           ! if physics grid resolution is not equal to fvm resolution
           call dyn2phys_all_vars(ie,                                         &
                ! spectral element state
                elem(ie)%state%dp3d(:,:,:,tl_f),                              &
                elem(ie)%state%T(:,:,:,tl_f),                                 &
                elem(ie)%derived%omega(:,:,:),                                &
                ! fvm state
                fvm(ie)%dp_fvm(:,:,:,n0_fvm),                                 &
                fvm(ie)%c(:,:,:,1:ntrac,n0_fvm),                              &
                pcnst, elem(ie)%metdet, fvm(ie),                              &
                !
                hyai(1)*ps0,                                                  &
                ! output
                dp3d_tmp(:,:,ie), ps_tmp(:,ie), q_tmp(:,:,:,ie),              &
                T_tmp(:,:,ie), omega_tmp(:,:,ie), phis_tmp(:,ie)              &
                )
           uv_tmp(:,:,:,ie) = &
                dyn2phys_vector(elem(ie)%state%v(:,:,:,:,tl_f),elem(ie))
         end do
       else
        ps_tmp(:,:)      = 0._r8
        T_tmp(:,:,:)     = 0._r8
        uv_tmp(:,:,:,:)  = 0._r8
        omega_tmp(:,:,:) = 0._r8
        phis_tmp(:,:)    = 0._r8
        Q_tmp(:,:,:,:)   = 0._r8
     endif ! iam < par%nprocs
     call t_stopf('dyn2phys')

     hybrid = config_thread_region(par,'serial')
     call get_loop_ranges(hybrid,ibeg=nets,iend=nete)

     ! Grab halo of surface fields (assuming we will use them)
     allocate(fld_phys_2d(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,1,halo_num_surface,nets:nete))
     call t_startf('dpcopy')
     do ie = nets, nete
        do j = 1, nc
           do i = 1, nc
              ioff = ((j-1)*nc) + i
              fld_phys_2d(i,j,1,halo_index_PHIS,ie) = phis_tmp(ioff, ie)
           end do
        end do
     end do
     call t_stopf('dpcopy')
     call Q3D_fill_halo(hybrid, elem, fld_phys_2d, nets, nete, 1, halo_num_surface, fvm)
     call t_startf('dpcopy')
     do index = 1, size(elem_fill_index)
        ie = elem_fill_index(index)%ie
        i = elem_fill_index(index)%ioff
        j = elem_fill_index(index)%joff
        phys_outs(halo_index_PHIS,  index) = fld_phys_2d(i,j,1,halo_index_PHIS,ie)
     end do
     do index = 1, size(ghost_fill_index)
        ie = ghost_fill_index(index)%ie
        i = ghost_fill_index(index)%ioff
        j = ghost_fill_index(index)%joff
        ioff = index + size(elem_fill_index)
        phys_outs(halo_index_PHIS, ioff) = fld_phys_2d(i,j,1,halo_index_PHIS,ie)
     end do
     call t_stopf('dpcopy')


     ! Grab halo of midpoint fields
     allocate(fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev,halo_num_midpoint,nets:nete))
     call t_startf('dpcopy')
     do ie = nets, nete
        do k = 1, nlev
           do j = 1, nc
              do i = 1, nc
                 ioff = ((j-1)*nc) + i
                 fld_phys(i,j,k,halo_index_u_con,ie) = uv_tmp(ioff,1,k,ie)
                 fld_phys(i,j,k,halo_index_v_con,ie) = uv_tmp(ioff,2,k,ie)
                 fld_phys(i,j,k,halo_index_Temp,ie)  = T_tmp(ioff,k,ie)
                 fld_phys(i,j,k,halo_index_OMEGA,ie) = omega_tmp(ioff,k,ie)
                 fld_phys(i,j,k,halo_index_DPDRY,ie) = dp3d_tmp(ioff,k,ie)
                 do m = 1, pcnst
                    fld_phys(i,j,k,halo_index_Q+m-1,ie) = q_tmp(ioff,k,m,ie)
                 end do
              end do
           end do
        end do
      end do
      call t_stopf('dpcopy')
!      do ie = nets, nete
!       call check(fld_phys(1:nc,1:nc,:,halo_index_u_con,ie),fv_nphys,0,0,nlev,-600.0_r8,600.0_r8,"beforeU")!DEBUGGING
!       call check(fld_phys(1:nc,1:nc,:,halo_index_v_con,ie),fv_nphys,0,0,nlev,-600.0_r8,600.0_r8,"beforeV")!DEBUGGING
!       call check(fld_phys(1:nc,1:nc,:,halo_index_Temp ,ie),fv_nphys,0,0,nlev,50.0_r8,600.0_r8,"beforeT")!DEBUGGING
!       call check(fld_phys(1:nc,1:nc,:,halo_index_OMEGA,ie),fv_nphys,0,0,nlev,-100.0_r8,100.0_r8,"beforeO")!DEBUGGING
!       do m = 1, pcnst
!         call check(fld_phys(1:nc,1:nc,:,halo_index_Q+m-1,ie),fv_nphys,0,0,nlev,-1.0E-12_r8,1.0_r8,"beforeQ")!DEBUGGING
!       end do
!     end do
     ! NB: istart_vector is location of U in fourth dim of fld_phys
     call Q3D_fill_halo(hybrid, elem, fld_phys, nets, nete, nlev, halo_num_midpoint, fvm, &
          istart_vector=halo_index_u_con, positive_definite_in=positive_definite)

     call t_startf('dpcopy')

!     do ie = nets, nete
!       call check(fld_phys(:,:,:,halo_index_u_con,ie),fv_nphys,nhc_phys,nhr_phys,nlev,-600.0_r8,600.0_r8,"U")!DEBUGGING
!       call check(fld_phys(:,:,:,halo_index_v_con,ie),fv_nphys,nhc_phys,nhr_phys,nlev,-600.0_r8,600.0_r8,"V")!DEBUGGING
!       call check(fld_phys(:,:,:,halo_index_Temp ,ie),fv_nphys,nhc_phys,nhr_phys,nlev,50.0_r8,600.0_r8,"T")!DEBUGGING
!       call check(fld_phys(:,:,:,halo_index_OMEGA,ie),fv_nphys,nhc_phys,nhr_phys,nlev,-100.0_r8,100.0_r8,"O")!DEBUGGING
!       do m = 1, pcnst
!         call check(fld_phys(:,:,:,halo_index_Q+m-1,ie),fv_nphys,nhc_phys,nhr_phys,nlev,-1.0E-12_r8,1.0_r8,"Q")!DEBUGGING
!       end do
!     end do
     do index = 1, size(elem_fill_index)
        ie = elem_fill_index(index)%ie
        i = elem_fill_index(index)%ioff
        j = elem_fill_index(index)%joff
        do k = 1, nlev
           phys_outm(k,halo_index_u_con,index) = fld_phys(i,j,k,halo_index_u_con,ie)
           phys_outm(k,halo_index_v_con,index) = fld_phys(i,j,k,halo_index_v_con,ie)
           phys_outm(k,halo_index_Temp, index) = fld_phys(i,j,k,halo_index_Temp,ie)
           phys_outm(k,halo_index_OMEGA,index) = fld_phys(i,j,k,halo_index_OMEGA,ie)
           do m = 1, pcnst
              phys_outm(k,halo_index_Q+m-1,index) = fld_phys(i,j,k,halo_index_Q+m-1,ie)
           end do
        end do
     end do
     do index = 1, size(ghost_fill_index)
        ie = ghost_fill_index(index)%ie
        i = ghost_fill_index(index)%ioff
        j = ghost_fill_index(index)%joff
        ioff = index + size(elem_fill_index)
        do k = 1, nlev
           phys_outm(k,halo_index_u_con,ioff) = fld_phys(i,j,k,halo_index_u_con,ie)
           phys_outm(k,halo_index_v_con,ioff) = fld_phys(i,j,k,halo_index_v_con,ie)
           phys_outm(k,halo_index_Temp, ioff) = fld_phys(i,j,k,halo_index_Temp,ie)
           phys_outm(k,halo_index_OMEGA,ioff) = fld_phys(i,j,k,halo_index_OMEGA,ie)
           do m = 1, pcnst
              phys_outm(k,halo_index_Q+m-1,ioff) = fld_phys(i,j,k,halo_index_Q+m-1,ie)
           end do
        end do
     end do
     call t_stopf('dpcopy')

     ! Grab halo of interface fields
     allocate(fld_phys_int(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev+1,halo_num_interface,nets:nete))
     do ie = nets, nete
       !
       ! compute full pressure (only including water vapor) at interfaces
       !
       fld_phys_int(:,:,1,halo_index_pint,ie) = hyai(1) * ps0
       do k = 1, nlev
         do j = 1, nc
           do i = 1, nc
             fld_phys_int(i,j,k+1,halo_index_pint,ie) = fld_phys_int(i,j,k,halo_index_pint,ie)+ &
                  (1.0_r8+fld_phys(i,j,k,halo_index_Q,ie))*fld_phys(i,j,k,halo_index_DPDRY,ie)
           end do
         end do
       end do
       call get_z_interface(&
            fld_phys_int(1:nc,1:nc,:,halo_index_pint,ie),&!dry pressure at interfaces
            fld_phys    (1:nc,1:nc,:,halo_index_Temp,ie),&!temperature (layer mean)
            fld_phys    (1:nc,1:nc,:,halo_index_Q   ,ie),&!dry water vapor mixing ratio (layer mean)
            fld_phys_2d (1:nc,1:nc,1,halo_index_PHIS,ie),&!surface geopotential
            fld_phys_int(1:nc,1:nc,:,halo_index_zi  ,ie)) !height at interfaces
     end do
     deallocate(fld_phys)
     deallocate(fld_phys_2d)
     call Q3D_fill_halo(hybrid, elem, fld_phys_int, nets, nete, nlev+1, halo_num_interface, fvm)
     call t_startf('dpcopy')

!     do ie = nets, nete
!       call check(fld_phys_int(:,:,:,halo_index_pint,ie),fv_nphys,nhc_phys,nhr_phys,nlev+1,0.0001_r8,200000.0_r8,"P")!DEBUGGING
!       call check(fld_phys_int(:,:,:,halo_index_zi  ,ie),fv_nphys,nhc_phys,nhr_phys,nlev+1,-1000.0_r8,100000.0_r8,"Z")!DEBUGGING
!     end do
     
     do index = 1, size(elem_fill_index)
        ie = elem_fill_index(index)%ie
        i = elem_fill_index(index)%ioff
        j = elem_fill_index(index)%joff
        do k = 1, nlev+1
          phys_outi(k, halo_index_pint,index) = fld_phys_int(i,j,k,halo_index_pint,ie)
          phys_outi(k, halo_index_zi  ,index) = fld_phys_int(i,j,k,halo_index_zi  ,ie)          
        end do
      end do
      do index = 1, size(ghost_fill_index)
        ie = ghost_fill_index(index)%ie
        i = ghost_fill_index(index)%ioff
        j = ghost_fill_index(index)%joff
        ioff = index + size(elem_fill_index)
        do k = 1, nlev+1
          phys_outi(k, halo_index_pint,ioff) = fld_phys_int(i,j,k,halo_index_pint,ie)
          phys_outi(k, halo_index_zi  ,ioff) = fld_phys_int(i,j,k,halo_index_zi  ,ie)          
        end do
      end do
     call t_stopf('dpcopy')
     deallocate(fld_phys_int)
     ! Deallocate the temporary arrays
     deallocate(ps_tmp)
     deallocate(dp3d_tmp)
     deallocate(phis_tmp)
     deallocate(T_tmp)
     deallocate(uv_tmp)
     deallocate(q_tmp)
     deallocate(omega_tmp)
   end subroutine q3d_halo_data

   subroutine q3d_halo_coord(elem_fill_index, ghost_fill_index, coord_out)
     use spmd_dyn,               only: local_dp_map
     use dimensions_mod,         only: nhc_phys
     use dyn_grid,               only: elem, fvm
     use coordinate_systems_mod, only: cart2spherical
     use coordinate_systems_mod, only: spherical_polar_t
     use shr_const_mod,          only: pi => shr_const_pi
     ! arguments
     type(element_index_t), intent(in)    :: elem_fill_index(:)
     type(halo_index_t),    intent(in)    :: ghost_fill_index(:)
     real(kind=r8),         intent(inout) :: coord_out(:,:)

 ! local variables
     integer                              :: nets, nete
     integer                              :: ie, i, j, index, ioff
     real(kind=r8),    allocatable        :: fld_phys_2d(:,:,:,:,:)
     type(hybrid_t)                       :: hybrid

     character(len=*), parameter          :: sub = 'q3d_halo_coord'

     real(kind=r8)                        :: dalpha, centerx, centery
     type (spherical_polar_t)             :: center_cart, sphere
!----------------------------------------------------------------------------

     if (fv_nphys /= 3) then
        call endrun(sub//': Q3D requires fv_nphys = 3')
     end if
     if (.not. local_dp_map) then
        call endrun(sub//': Q3D only supports local maps')
     end if

     hybrid = config_thread_region(par,'serial')
     call get_loop_ranges(hybrid,ibeg=nets,iend=nete)

     ! Grab halo of coordinates
     allocate(fld_phys_2d(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,1,2,nets:nete))
     call t_startf('dpcopy')
     do ie = nets, nete
       dalpha  = abs(elem(ie)%corners(1)%x-elem(ie)%corners(2)%x)/nc
       fld_phys_2d(:,:,:,:,ie) = HUGE(1.0_r8)
       do j = 1-nhc_phys, nc+nhc_phys
         do i = 1-nhc_phys, nc+nhc_phys
           centerx = tan(elem(ie)%corners(1)%x+(i-0.5_r8)*dalpha)
           centery = tan(elem(ie)%corners(1)%y+(j-0.5_r8)*dalpha)
           sphere  = cart2spherical(centerx,centery,elem(ie)%FaceNum)
           ioff = ((j-1)*nc) + i            !not used???
           fld_phys_2d(i,j,1,halo_index_lon,ie) = sphere%lon
           fld_phys_2d(i,j,1,halo_index_lat,ie) = sphere%lat
         end do
       end do
     end do
     !! PHL: Replace below with analytic off-face coordinates
!     call Q3D_fill_halo(hybrid, elem, fld_phys_2d, nets, nete, 1, 2, fvm)
     call t_startf('dpcopy')
     do index = 1, size(elem_fill_index)
        ie = elem_fill_index(index)%ie
        i = elem_fill_index(index)%ioff
        j = elem_fill_index(index)%joff
        coord_out(halo_index_lon, index) = fld_phys_2d(i,j,1,halo_index_lon,ie)
        coord_out(halo_index_lat, index) = fld_phys_2d(i,j,1,halo_index_lat,ie)
     end do
     do index = 1, size(ghost_fill_index)
        ie = ghost_fill_index(index)%ie
        i = ghost_fill_index(index)%ioff
        j = ghost_fill_index(index)%joff
        ioff = index + size(elem_fill_index)
        coord_out(halo_index_lon, ioff) = fld_phys_2d(i,j,1,halo_index_lon,ie)
        coord_out(halo_index_lat, ioff) = fld_phys_2d(i,j,1,halo_index_lat,ie)
     end do
     call t_stopf('dpcopy')
     deallocate(fld_phys_2d)
   end subroutine q3d_halo_coord

   subroutine q3d_set_indexing_ghost(ghost_fill_index)
     use dyn_grid,               only: elem

     type(halo_index_t),    intent(inout) :: ghost_fill_index(:)

     integer                              :: index, ie, i, j, k
     integer                              :: ighost, jghost
     integer                              :: gcid, ghi, ghj
     character(len=128)                   :: errmsg
     character(len=*), parameter          :: sub = 'q3d_set_indexing_ghost'

     do index = 1, size(ghost_fill_index)
        ghost_fill_index(index)%ie = 0
        do ie = 1, nelemd
           k = 1
           do j = 1, fv_nphys
              do i = 1, fv_nphys
                 gcid = ((elem(ie)%GlobalId - 1)*fv_nphys*fv_nphys) + k
                 if (gcid == ghost_fill_index(index)%gcid) then
                    ghost_fill_index(index)%ie = ie
                    ! Set ioff and joff by checking actual edges
                    if (elem(ie)%FaceNum > 4) then
                       ! Top and bottom are rotated
                       ighost = abs(ghost_fill_index(index)%jghost)
                       jghost = abs(ghost_fill_index(index)%ighost)
                    else
                       ighost = abs(ghost_fill_index(index)%ighost)
                       jghost = abs(ghost_fill_index(index)%jghost)
                    end if
                    if (i == 1) then
                       if (j == 1) then
                          ghi = 1 - ighost
                          ghj = 1 - jghost
                       else if (j == fv_nphys) then
                          ghi = 1 - ighost
                          ghj = fv_nphys + jghost
                       else
                          ghi = 1 - ighost
                          if (jghost /= 0) then
                             write(errmsg, '(a,4(a,i0))') ': Bad jghost with i = 1', &
                                  ', j = ',j,', face = ',elem(ie)%FaceNum, &
                                  ', ghost = ',ighost,', ',jghost
                             call endrun(sub//trim(errmsg))
                          end if
                          ghj = j
                       end if
                    else if (i == fv_nphys) then
                       if (j == 1) then
                          ghi = fv_nphys + ighost
                          ghj = 1 - jghost
                       else if (j == fv_nphys) then
                          ghi = fv_nphys + ighost
                          ghj = fv_nphys + jghost
                       else
                          ghi = fv_nphys + ighost
                          if (jghost /= 0) then
                             write(errmsg, '(a,4(a,i0))') ': Bad jghost with i = fv_nphys', &
                                  ', j = ',j,', face = ',elem(ie)%FaceNum, &
                                  ', ghost = ',ighost,', ',jghost
                             call endrun(sub//trim(errmsg))
                          end if
                          ghj = j
                       end if
                    else if (j == 1) then
                       ghi = i
                       if (ighost /= 0) then
                          write(errmsg, '(a,4(a,i0))') ': Bad ighost with j = 1', &
                               ', i = ',i,', face = ',elem(ie)%FaceNum,      &
                               ', ghost = ',ighost,', ',jghost
                          call endrun(sub//trim(errmsg))
                       end if
                       ghj = 1 - jghost
                    else if (j == fv_nphys) then
                       ghi = i
                       if (ighost /= 0) then
                          write(errmsg, '(a,4(a,i0))') ': Bad ighost with j = fv_nphys', &
                               ', i = ',i,', face = ',elem(ie)%FaceNum,       &
                               ', ghost = ',ighost,', ',jghost
                          call endrun(sub//trim(errmsg))
                       end if
                       ghj = fv_nphys + jghost
                    else
                       write(errmsg, '(a,i0,a,i0,a)') ': ERROR: Interior point (',i,', ',j,')'
                       call endrun(sub//trim(errmsg))
                    end if
                    ghost_fill_index(index)%ioff = ghi
                    ghost_fill_index(index)%joff = ghj
                    exit
                 end if
                 k = k + 1
              end do
              if (ghost_fill_index(index)%ie > 0) then
                 exit
              end if
           end do
           if (ghost_fill_index(index)%ie > 0) then
              exit
           end if
        end do
        if (ghost_fill_index(index)%ie < 1) then
           write(errmsg, '(a,i0)') ': global column not found for gcid = ',ghost_fill_index(index)%gcid
           call endrun(sub//trim(errmsg))
        end if
     end do

   end subroutine q3d_set_indexing_ghost

   subroutine q3d_set_indexing_element(elem_fill_index)
     use dyn_grid,               only: elem

     type(element_index_t), intent(inout) :: elem_fill_index(:)

     integer                              :: index, ie, i, j, k
     integer                              :: gcid, ghi, ghj
     character(len=128)                   :: errmsg
     character(len=*), parameter          :: sub = 'q3d_set_indexing_element'

     do index = 1, size(elem_fill_index)
        elem_fill_index(index)%ie = 0
        do ie = 1, nelemd
           k = 1
           do j = 1, fv_nphys
              do i = 1, fv_nphys
                 gcid = ((elem(ie)%GlobalId - 1)*fv_nphys*fv_nphys) + k
                 if (gcid == elem_fill_index(index)%gcid) then
                    elem_fill_index(index)%ie = ie
                    elem_fill_index(index)%ioff = i
                    elem_fill_index(index)%joff = j
                 end if
                 k = k + 1
              end do
              if (elem_fill_index(index)%ie > 0) then
                 exit
              end if
           end do
           if (elem_fill_index(index)%ie > 0) then
              exit
           end if
        end do
        if (elem_fill_index(index)%ie < 1) then
           write(errmsg, '(a,i0)') ': global column not found for gcid = ',elem_fill_index(index)%gcid
           call endrun(sub//trim(errmsg))
        end if
     end do

   end subroutine q3d_set_indexing_element

   subroutine Q3D_fill_halo(hybrid, elem, fld_phys, nets, nete, num_lev, num_flds, fvm, istart_vector, positive_definite_in)
     use element_mod,            only: element_t
     use fvm_control_volume_mod, only: fvm_struct
     use dimensions_mod,         only: np, nhc_phys, fv_nphys, nhr_phys, ns_phys
     use hybrid_mod,             only: hybrid_t
     use fvm_reconstruction_mod, only: extend_panel_interpolate
     use fvm_mapping,            only: fill_halo_phys

     type (hybrid_t), intent(in)   :: hybrid  ! distributed parallel structure (shared)
     integer       , intent(in)    :: nets,nete,num_flds,num_lev
     real (kind=r8), intent(inout) :: fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys, &
          num_lev,num_flds, nets:nete)
     type (element_t)     , intent(inout) :: elem(:)
     type(fvm_struct)     , intent(in)    :: fvm(:)
     integer, optional    , intent(in)    :: istart_vector
     logical, optional    , intent(in)    :: positive_definite_in(num_flds)

     integer          :: i, j, ie, k, itr
     real (kind=r8)   :: v(2)
     real (kind=r8)   :: fld_phys_out(1-nhr_phys:fv_nphys+nhr_phys,1-nhr_phys:fv_nphys+nhr_phys)
     logical          :: positive_definite(num_flds)
     real (kind=r8)   :: D(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,2,2,nets:nete)
     real (kind=r8)   :: fill_value = -1.0_r8

     call get_extend_panel_contra_to_latlon_matrix(elem, fvm, hybrid, nets, nete,D)
     
     if (present(positive_definite_in)) then
        positive_definite = positive_definite_in
     else
        positive_definite = .false.
     end if

     call fill_halo_phys(elem,fld_phys,hybrid,nets,nete,num_lev,num_flds)

     if (present(istart_vector)) then
       do ie = nets, nete
         do k = 1, num_lev
           do j = 1-nhr_phys, fv_nphys+nhr_phys
             do i = 1-nhr_phys, fv_nphys+nhr_phys
               !
               ! convert lat-lon vectors to contra-variant gnomonic
               !
               v(1) = fld_phys(i,j,k,istart_vector  ,ie)
               v(2) = fld_phys(i,j,k,istart_vector+1,ie)
               fld_phys(i,j,k,istart_vector  ,ie)=fvm(ie)%Dinv_physgrid(i,j,1,1)*v(1) + fvm(ie)%Dinv_physgrid(i,j,1,2)*v(2)
               fld_phys(i,j,k,istart_vector+1,ie)=fvm(ie)%Dinv_physgrid(i,j,2,1)*v(1) + fvm(ie)%Dinv_physgrid(i,j,2,2)*v(2)
             end do
           end do
         end do
       end do
     end if
     !
     ! third row halo not used as well as corners in 2 deep halo - set values to 9E99
     !
     do ie = nets, nete
       do itr = 1, num_flds
         do k = 1, num_lev
           call extend_panel_interpolate(fv_nphys,nhc_phys,nhr_phys,nhr_phys,ns_phys,nhr_phys,&
                fld_phys(:,:,k,itr,ie),fvm(ie)%cubeboundary,fvm(ie)%halo_interp_weight,fvm(ie)%ibase,&
                fld_phys_out(:,:))
           !
           ! set unfilled corner values
           !
           fld_phys_out(1-nhr_phys       ,1-nhr_phys       ) = fill_value
           fld_phys_out(fv_nphys+nhr_phys,1-nhr_phys       ) = fill_value
           fld_phys_out(fv_nphys+nhr_phys,fv_nphys+nhr_phys) = fill_value
           fld_phys_out(1-nhr_phys       ,fv_nphys+nhr_phys) = fill_value
           !
           ! set third row halo (not used by channels)
           !
           fld_phys(1-nhc_phys       ,:,k,itr,ie) = fill_value
           fld_phys(fv_nphys+nhc_phys,:,k,itr,ie) = fill_value
           fld_phys(:,1-nhc_phys       ,k,itr,ie) = fill_value
           fld_phys(:,fv_nphys+nhc_phys,k,itr,ie) = fill_value           
           fld_phys(1-nhr_phys:fv_nphys+nhr_phys,1-nhr_phys:fv_nphys+nhr_phys,k,itr,ie) = fld_phys_out(:,:)
         end do
         if (positive_definite(itr)) then
           do k = 1, num_lev
             do j = 1-nhr_phys, fv_nphys+nhr_phys
               do i = 1-nhr_phys, fv_nphys+nhr_phys
                 fld_phys(i,j,k,itr,ie) = MAX(fld_phys(i,j,k,itr,ie),0.0_r8)
               end do
             end do
           end do
         end if
       end do
     end do

     if (present(istart_vector)) then
        do ie = nets, nete
           do k = 1, num_lev
              do j = 1-nhr_phys, fv_nphys+nhr_phys
                 do i = 1-nhr_phys, fv_nphys+nhr_phys
                    !
                    ! convert lat-lon vectors to contra-variant gnomonic
                    !
                    v(1) = fld_phys(i,j,k,istart_vector  ,ie)
                    v(2) = fld_phys(i,j,k,istart_vector+1,ie)
                    fld_phys(i,j,k,istart_vector  ,ie)=D(i,j,1,1,ie)*v(1) + D(i,j,1,2,ie)*v(2)
                    fld_phys(i,j,k,istart_vector+1,ie)=D(i,j,2,1,ie)*v(1) + D(i,j,2,2,ie)*v(2)
                 end do
              end do
           end do
        end do
     end if

      
   end subroutine Q3D_fill_halo

   subroutine get_extend_panel_contra_to_latlon_matrix(elem, fvm, hybrid, nets, nete,D)
     use element_mod,            only: element_t
     use fvm_control_volume_mod, only: fvm_struct
     use dimensions_mod,         only: fv_nphys, nhc_phys
     use cube_mod               ,only: dmap
     use control_mod,            only: cubed_sphere_map
     
     type (element_t) , intent(in)    :: elem(:)
     type (fvm_struct), intent(in)    :: fvm(:)
     type (hybrid_t)  , intent(in)    :: hybrid
     
     integer, intent(in)              :: nets,nete
     ! D is derivative of gnomonic mapping
     real (kind=r8), intent(out)      :: D(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,2,2,nets:nete)
     
     integer                          :: ie, i, j     
     real (kind=r8)                   :: x(2)
     !
     ! create a normalized element coordinate system with a halo
     !    
     do ie=nets,nete
       do j=1-nhc_phys,fv_nphys+nhc_phys
         do i=1-nhc_phys,fv_nphys+nhc_phys          
           x(1) = elem(ie)%corners(1)%x+(i-0.5_r8)*fvm(ie)%dalpha_physgrid
           x(2) = elem(ie)%corners(1)%y+(j-0.5_r8)*fvm(ie)%dalpha_physgrid
           !
           ! convert to element normalized coordinates
           !
           x(1) =(x(1)-elem(ie)%corners(1)%x)/&
                (0.5_r8*dble(fv_nphys)*fvm(ie)%dalpha_physgrid)-1.0_r8
           x(2) =(x(2)-elem(ie)%corners(1)%y)/&
                (0.5_r8*dble(fv_nphys)*fvm(ie)%dalpha_physgrid)-1.0_r8
           !
           ! compute D
           !
           call Dmap(D(i,j,:,:,ie),x(1),x(2),elem(ie)%corners3D,cubed_sphere_map,elem(ie)%corners,elem(ie)%u2qmap,elem(ie)%facenum)
         end do
       end do
     end do
   end subroutine get_extend_panel_contra_to_latlon_matrix
   
  subroutine get_z_interface(pint,t,qwater,phis,zint)
    use dimensions_mod, only: nlev,nc
    use physconst,      only: rair, epsilo, gravit
    use shr_vmath_mod,  only: shr_vmath_log
    implicit none
    
    real(r8), intent(in)  :: pint  (nc,nc,nlev+1)  ! full pressure at interfaces
    real(r8), intent(in)  :: qwater(nc,nc,nlev)    ! dry water vapor mixing ratio
    real(r8), intent(in)  :: t     (nc,nc,nlev)    ! temperature
    real(r8), intent(in)  :: phis  (nc,nc)         ! surface geopotential
    
    real(r8), intent(out) :: zint(nc,nc,nlev+1)    ! height at interface
    !
    !---------------------------Local variables-----------------------------
    !
    integer  :: i,j,k
    real(r8) :: hkl(nc,nc)            ! off-diagonal element
    real(r8) :: rog                   ! Rair / gravit
    real(r8) :: tv(nc,nc,nlev)        ! Virtual temperature
    real(r8) :: log_pint(nc,nc,nlev+1)! log of pint
    real(r8) :: inv_epsilon
    
    rog = rair/gravit
    inv_epsilon = 1/Epsilo
    !
    ! compute virtual temperature
    !
    do k=1,nlev
      do j=1,nc
        do i=1,nc
          tv(i,j,k)  = T(i,j,k)*(1.0_r8+inv_epsilon*qwater(i,j,k))/(1.0_r8+qwater(i,j,k))
        end do
      end do
    end do
    !
    ! compute log of pint
    !
    call shr_vmath_log(pint,log_pint,nc*nc*(nlev+1)) 
    !
    ! Compute zi from bottom up.
    !
    zint(:,:,nlev+1) = phis(:,:)/gravit
    do k = nlev, 1, -1
      hkl(:,:) = (log_pint(:,:,k+1) - log_pint(:,:,k))
      ! dynamics uses hkl(:,:) = pdel(:,:,k) / pmid(:,:,k)
      zint(:,:,k) = zint(:,:,k+1) + rog * tv(:,:,k) * hkl(:,:)
    end do
  end subroutine get_z_interface


  subroutine check(f,nc,nhc,nhr,nlev,min,max,str)
    implicit none
    integer, intent(in)   :: nlev,nc,nhr,nhc    
    real(r8), intent(in)  :: f(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev)  ! full pressure at interfaces
    real(r8), intent(in)  :: min,max
    character(len=*):: str
    integer:: i,j,k
    do k=1,nlev
      !
      ! DO NOT CHECK THIRD ROW HALO
      !
      if (nhr.ne.0) then
        do j=1-nhr,nc+nhr
          do i=1-nhr,nc+nhr
            !
            ! DO NOT CHECK CORNERS 
            !
            if ((i.ne.1-nhr .and.(j.ne.1-nhr.or.j.ne.nc+nhr)).and.&
                 (i.ne.nc+nhr.and.(j.ne.1-nhr.or.j.ne.nc+nhr))) then
              if (f(i,j,k)<min.or.f(i,j,k)>max) then
                write(*,*) "xxx ERROR",str,i,j,k,f(i,j,k),min,max
              end if
            end if
          end do
        end do
      else
        do j=1,nc
          do i=1,nc
            if (f(i,j,k)<min.or.f(i,j,k)>max) then
              write(*,*) "xxx ERROR",str,i,j,k,f(i,j,k),min,max
            end if
          end do
        end do
      end if
    end do    
  end subroutine check
end module q3d_coupling
