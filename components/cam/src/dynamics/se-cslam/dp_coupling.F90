module dp_coupling

!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
use constituents,   only: pcnst, cnst_type

use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
use spmd_utils,     only: mpicom, iam
use dyn_grid,       only: get_gcol_block_d, TimeLevel, edgebuf
use dyn_comp,       only: dyn_export_t, dyn_import_t

use physics_types,  only: physics_state, physics_tend
use phys_grid,      only: get_ncols_p, get_gcol_all_p, block_to_chunk_send_pters, &
                          transpose_block_to_chunk, block_to_chunk_recv_pters,    &
                          chunk_to_block_send_pters, transpose_chunk_to_block,    &
                          chunk_to_block_recv_pters
use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk, pbuf_get_field

use dp_mapping,     only: nphys_pts

use cam_logfile,    only: iulog
use perf_mod,       only: t_startf, t_stopf, t_barrierf
use cam_abortutils, only: endrun

use parallel_mod,   only: par
use thread_mod,     only: horz_num_threads
use hybrid_mod,     only: config_thread_region, get_loop_ranges, hybrid_t
use dimensions_mod, only: np, npsq, nelemd, nlev, nc, qsize, ntrac, fv_nphys

use dof_mod,        only: UniquePoints, PutUniquePoints
use element_mod,    only: element_t
use fvm_control_volume_mod, only: fvm_struct

implicit none
private
save

public :: d_p_coupling, p_d_coupling

real (kind=r8), allocatable :: q_prev(:,:,:,:) ! Previous Q for computing tendencies

!=========================================================================================
CONTAINS
!=========================================================================================

subroutine d_p_coupling(phys_state, phys_tend,  pbuf2d, dyn_out)

   use gravity_waves_sources,  only: gws_src_fnct
   use dyn_comp,               only: frontgf_idx, frontga_idx
   use phys_control,           only: use_gw_front, use_gw_front_igw
   use hycoef,                 only: hyai, ps0
   use fvm_control_volume_mod, only: n0_fvm
   use fvm_mapping,            only: dyn2phys_vector, dyn2phys_all_vars
   use time_mod,               only: timelevel_qdp
   use control_mod,            only: qsplit
   use test_fvm_mapping,       only: test_mapping_overwrite_dyn_state, test_mapping_output_phys_state

   ! arguments
   type(dyn_export_t),  intent(inout)                               :: dyn_out             ! dynamics export
   type(physics_buffer_desc), pointer                               :: pbuf2d(:,:)
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend


   ! LOCAL VARIABLES
   type(element_t), pointer     :: elem(:)             ! pointer to dyn_out element array
   integer                      :: ie                  ! indices over elements
   integer                      :: lchnk, icol, ilyr   ! indices over chunks, columns, layers

   real (kind=r8),  allocatable :: ps_tmp(:,:)         ! temp array to hold ps
   real (kind=r8),  allocatable :: dp3d_tmp(:,:,:)     ! temp array to hold dp3d
   real (kind=r8),  allocatable :: dp3d_tmp_tmp(:,:)
   real (kind=r8),  allocatable :: phis_tmp(:,:)       ! temp array to hold phis
   real (kind=r8),  allocatable :: T_tmp(:,:,:)        ! temp array to hold T
   real (kind=r8),  allocatable :: uv_tmp(:,:,:,:)     ! temp array to hold u and v
   real (kind=r8),  allocatable :: q_tmp(:,:,:,:)      ! temp to hold advected constituents
   real (kind=r8),  allocatable :: omega_tmp(:,:,:)    ! temp array to hold omega

   ! Frontogenesis
   real (kind=r8),  allocatable :: frontgf(:,:,:)      ! temp arrays to hold frontogenesis
   real (kind=r8),  allocatable :: frontga(:,:,:)      ! function (frontgf) and angle (frontga)
                                                        ! Pointers to pbuf
   real (kind=r8),  pointer     :: pbuf_frontgf(:,:)
   real (kind=r8),  pointer     :: pbuf_frontga(:,:)

   integer                      :: ncols,i,j,ierr,k,iv
   integer                      :: ioff, m, m_cnst
   integer                      :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)
   integer                      :: tsize               ! amount of data per grid point passed to physics
   integer,         allocatable :: bpter(:,:)          ! offsets into block buffer for packing data
   integer                      :: cpter(pcols,0:pver) ! offsets into chunk buffer for unpacking data
   integer                      :: nphys

   real (kind=r8),  allocatable :: bbuffer(:), cbuffer(:) ! transpose buffers
   real (kind=r8),  allocatable :: qgll(:,:,:,:)
   real (kind=r8)               :: inv_dp3d(np,np,nlev)
   integer                      :: tl_f, tl_qdp_np0, tl_qdp_np1
   logical                      :: lmono

   type(physics_buffer_desc), pointer :: pbuf_chnk(:)
   !----------------------------------------------------------------------------

   elem => dyn_out%elem
   tl_f = TimeLevel%n0
   call TimeLevel_Qdp(TimeLevel, qsplit, tl_qdp_np0,tl_qdp_np1)

   nullify(pbuf_chnk)
   nullify(pbuf_frontgf)
   nullify(pbuf_frontga)

   if (fv_nphys > 0) then
      nphys = fv_nphys
   else
     allocate(qgll(np,np,nlev,pcnst))
     nphys = np
   end if

   ! Allocate temporary arrays to hold data for physics decomposition
   allocate(ps_tmp(nphys_pts,nelemd))
   allocate(dp3d_tmp(nphys_pts,pver,nelemd))
   allocate(dp3d_tmp_tmp(nphys_pts,pver))
   allocate(phis_tmp(nphys_pts,nelemd))
   allocate(T_tmp(nphys_pts,pver,nelemd))
   allocate(uv_tmp(nphys_pts,2,pver,nelemd))
   allocate(q_tmp(nphys_pts,pver,pcnst,nelemd))
   allocate(omega_tmp(nphys_pts,pver,nelemd))

   if (use_gw_front .or. use_gw_front_igw) then
      allocate(frontgf(nphys_pts,pver,nelemd), stat=ierr)
      if (ierr /= 0) call endrun("dp_coupling: Allocate of frontgf failed.")
      allocate(frontga(nphys_pts,pver,nelemd), stat=ierr)
      if (ierr /= 0) call endrun("dp_coupling: Allocate of frontga failed.")
   end if

   if (iam < par%nprocs) then
      if (use_gw_front .or. use_gw_front_igw) then
         call gws_src_fnct(elem, tl_f, tl_qdp_np0, frontgf, frontga, nphys)
      end if

      if (fv_nphys > 0) then
         call test_mapping_overwrite_dyn_state(elem,dyn_out%fvm)
         !******************************************************************
         ! physics runs on an FVM grid: map GLL vars to physics grid
         !******************************************************************
         call t_startf('dyn2phys')
         do ie = 1, nelemd
            ! note that the fvm halo has been filled in prim_run_subcycle
            ! if physics grid resolution is not equal to fvm resolution
            call dyn2phys_all_vars(ie,                                         &
               ! spectral element state
               elem(ie)%state%dp3d(:,:,:,tl_f),                                &
               elem(ie)%state%T(:,:,:,tl_f),                                   &
               elem(ie)%derived%omega(:,:,:),                                  &
               ! fvm state
               dyn_out%fvm(ie)%dp_fvm(:,:,:,n0_fvm),                           &
               dyn_out%fvm(ie)%c(:,:,:,1:ntrac,n0_fvm),                        &
               pcnst, elem(ie)%metdet, dyn_out%fvm(ie),        &
               !
               hyai(1)*ps0,                                                    &
               ! output
               dp3d_tmp(:,:,ie), ps_tmp(:,ie), q_tmp(:,:,:,ie), T_tmp(:,:,ie), &
               omega_tmp(:,:,ie), phis_tmp(:,ie)                               &
               )
            uv_tmp(:,:,:,ie) = &
               dyn2phys_vector(elem(ie)%state%v(:,:,:,:,tl_f),elem(ie))
         end do
         call t_stopf('dyn2phys')
      else

         !******************************************************************
         ! Physics runs on GLL grid: collect unique points before mapping to
         ! physics decomposition
         !******************************************************************

         if (qsize < pcnst) then
            call endrun('d_p_coupling: Fewer GLL tracers advected than required')
         end if

         call t_startf('UniquePoints')
         do ie = 1, nelemd
           inv_dp3d(:,:,:) = 1.0_r8/elem(ie)%state%dp3d(:,:,:,tl_f)
           do m=1,pcnst
             qgll(:,:,:,m) = elem(ie)%state%Qdp(:,:,:,m,tl_qdp_np0)*inv_dp3d(:,:,:)
           end do
            ncols = elem(ie)%idxP%NumUniquePts
            call UniquePoints(elem(ie)%idxP, elem(ie)%state%psdry(:,:), ps_tmp(1:ncols,ie))
            call UniquePoints(elem(ie)%idxP, nlev, elem(ie)%state%dp3d(:,:,:,tl_f), dp3d_tmp(1:ncols,:,ie))
            call UniquePoints(elem(ie)%idxP, nlev, elem(ie)%state%T(:,:,:,tl_f), T_tmp(1:ncols,:,ie))
            call UniquePoints(elem(ie)%idxV, 2, nlev, elem(ie)%state%V(:,:,:,:,tl_f), uv_tmp(1:ncols,:,:,ie))
            call UniquePoints(elem(ie)%idxV, nlev, elem(ie)%derived%omega, omega_tmp(1:ncols,:,ie))

            call UniquePoints(elem(ie)%idxP, elem(ie)%state%phis, phis_tmp(1:ncols,ie))
            call UniquePoints(elem(ie)%idxP, nlev, pcnst, qgll,Q_tmp(1:ncols,:,:,ie))
         end do
         call t_stopf('UniquePoints')

      end if ! if fv_nphys>0

   else

      ps_tmp(:,:)      = 0._r8
      T_tmp(:,:,:)     = 0._r8
      uv_tmp(:,:,:,:)  = 0._r8
      omega_tmp(:,:,:) = 0._r8
      phis_tmp(:,:)    = 0._r8
      Q_tmp(:,:,:,:)   = 0._r8

      if (use_gw_front .or. use_gw_front_igw) then
         frontgf(:,:,:) = 0._r8
         frontga(:,:,:) = 0._r8
      end if

   endif ! iam < par%nprocs

   if (fv_nphys<1) deallocate(qgll)

   ! q_prev is for saving the tracer fields for calculating tendencies
   if (.not. allocated(q_prev)) then
      allocate(q_prev(pcols,pver,pcnst,begchunk:endchunk))
   end if
   q_prev = 0.0_R8

   call t_startf('dpcopy')
   if (local_dp_map) then

      !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, ilyr, m, pbuf_chnk, pbuf_frontgf, pbuf_frontga)
      do lchnk = begchunk, endchunk

         ncols = get_ncols_p(lchnk)
         call get_gcol_all_p(lchnk, pcols, pgcols)

         pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
         if (use_gw_front .or. use_gw_front_igw) then
            call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
            call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
         end if

         do icol = 1, ncols
            call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
            ie   = idmb3(1)
            ioff = idmb2(1)
            phys_state(lchnk)%ps(icol)   = ps_tmp(ioff,ie)
            phys_state(lchnk)%phis(icol) = phis_tmp(ioff,ie)
            do ilyr=1,pver
               phys_state(lchnk)%pdel(icol,ilyr)  = dp3d_tmp(ioff,ilyr,ie)
               phys_state(lchnk)%t(icol,ilyr)     = T_tmp(ioff,ilyr,ie)
               phys_state(lchnk)%u(icol,ilyr)     = uv_tmp(ioff,1,ilyr,ie)
               phys_state(lchnk)%v(icol,ilyr)     = uv_tmp(ioff,2,ilyr,ie)
               phys_state(lchnk)%omega(icol,ilyr) = omega_tmp(ioff,ilyr,ie)

               if (use_gw_front .or. use_gw_front_igw) then
                  pbuf_frontgf(icol,ilyr) = frontgf(ioff,ilyr,ie)
                  pbuf_frontga(icol,ilyr) = frontga(ioff,ilyr,ie)
               endif
            end do

            do m = 1, pcnst
               do ilyr = 1, pver
                  phys_state(lchnk)%q(icol,ilyr,m) = Q_tmp(ioff,ilyr,m,ie)
               end do
            end do
         end do
      end do

   else  ! .not. local_dp_map

      tsize = 5 + pcnst
      if (use_gw_front .or. use_gw_front_igw) tsize = tsize + 2

      allocate(bbuffer(tsize*block_buf_nrecs))
      allocate(cbuffer(tsize*chunk_buf_nrecs))
      if (fv_nphys > 0) then
         allocate(bpter(fv_nphys*fv_nphys,0:pver))
      else
         allocate(bpter(npsq,0:pver))
      end if

      if (iam < par%nprocs) then
         !$omp parallel do num_threads(horz_num_threads) private (ie, bpter, icol, ilyr, m, ncols, ioff)
         do ie = 1, nelemd

            if (fv_nphys > 0) then
               call block_to_chunk_send_pters(elem(ie)%GlobalID, fv_nphys*fv_nphys,  &
                  pver+1, tsize, bpter)
               ncols = fv_nphys*fv_nphys
            else
               call block_to_chunk_send_pters(elem(ie)%GlobalID, npsq,       &
                  pver+1, tsize, bpter)
               ncols = elem(ie)%idxP%NumUniquePts
            end if

            do icol=1,ncols
               bbuffer(bpter(icol,0)+2:bpter(icol,0)+tsize-1) = 0.0_r8
               bbuffer(bpter(icol,0))   = ps_tmp(icol,ie)
               bbuffer(bpter(icol,0)+1) = phis_tmp(icol,ie)

               do ilyr=1,pver
                  ioff = 0
                  bbuffer(bpter(icol,ilyr)+ioff) = T_tmp(icol,ilyr,ie)
                  ioff = ioff + 1
                  bbuffer(bpter(icol,ilyr)+ioff) = uv_tmp(icol,1,ilyr,ie)
                  ioff = ioff + 1
                  bbuffer(bpter(icol,ilyr)+ioff) = uv_tmp(icol,2,ilyr,ie)
                  ioff = ioff + 1
                  bbuffer(bpter(icol,ilyr)+ioff) = omega_tmp(icol,ilyr,ie)
                  ioff = ioff + 1
                  bbuffer(bpter(icol,ilyr)+ioff) = dp3d_tmp(icol,ilyr,ie)
                  if (use_gw_front .or. use_gw_front_igw) then
                     ioff = ioff + 1
                     bbuffer(bpter(icol,ilyr)+ioff) = frontgf(icol,ilyr,ie)
                     ioff = ioff + 1
                     bbuffer(bpter(icol,ilyr)+ioff) = frontga(icol,ilyr,ie)
                  end if

                  do m=1,pcnst
                     bbuffer(bpter(icol,ilyr)+tsize-pcnst-1+m) = Q_tmp(icol,ilyr,m,ie)
                  end do
               end do
            end do
         end do

      else
         bbuffer(:) = 0._r8
      end if

      call t_barrierf ('sync_blk_to_chk', mpicom)
      call t_startf ('block_to_chunk')
      call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
      call t_stopf  ('block_to_chunk')

      !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, cpter, icol, ilyr, m, pbuf_chnk, pbuf_frontgf, pbuf_frontga, ioff)
      do lchnk = begchunk, endchunk
         ncols = phys_state(lchnk)%ncol

         pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

         if (use_gw_front .or. use_gw_front_igw) then
            call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
            call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
         end if

         call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

         do icol = 1, ncols
            phys_state(lchnk)%ps(icol)   = cbuffer(cpter(icol,0))
            phys_state(lchnk)%phis(icol) = cbuffer(cpter(icol,0)+1)

            do ilyr = 1, pver
               ioff = 0
               phys_state(lchnk)%t(icol,ilyr)     = cbuffer(cpter(icol,ilyr)+ioff)
               ioff = ioff + 1
               phys_state(lchnk)%u(icol,ilyr)     = cbuffer(cpter(icol,ilyr)+ioff)
               ioff = ioff + 1
               phys_state(lchnk)%v(icol,ilyr)     = cbuffer(cpter(icol,ilyr)+ioff)
               ioff = ioff + 1
               phys_state(lchnk)%omega(icol,ilyr) = cbuffer(cpter(icol,ilyr)+ioff)
               ioff = ioff + 1
               phys_state(lchnk)%pdel(icol,ilyr)  = cbuffer(cpter(icol,ilyr)+ioff)

               if (use_gw_front .or. use_gw_front_igw) then
                  ioff = ioff + 1
                  pbuf_frontgf(icol,ilyr) = cbuffer(cpter(icol,ilyr)+ioff)
                  ioff = ioff + 1
                  pbuf_frontga(icol,ilyr) = cbuffer(cpter(icol,ilyr)+ioff)
               endif

               do m = 1, pcnst
                  phys_state(lchnk)%q  (icol,ilyr,m) = cbuffer(cpter(icol,ilyr)+tsize-pcnst-1+m)
               end do

            end do
         end do
      end do

      deallocate( bbuffer )
      deallocate( cbuffer )

   end if
   call t_stopf('dpcopy')

   ! Save the tracer fields for calculating tendencies
   do lchnk = begchunk, endchunk
      ncols = phys_state(lchnk)%ncol
      q_prev(1:ncols,1:pver,1:pcnst,lchnk) = phys_state(lchnk)%q(1:ncols,1:pver,1:pcnst)
   end do
   call test_mapping_output_phys_state(phys_state,dyn_out%fvm)

   ! Deallocate the temporary arrays
   deallocate(ps_tmp)
   deallocate(dp3d_tmp)
   deallocate(phis_tmp)
   deallocate(T_tmp)
   deallocate(uv_tmp)
   deallocate(q_tmp)
   deallocate(omega_tmp)

   call t_startf('derived_phys')
   call derived_phys_dry(phys_state, phys_tend, pbuf2d)
   call t_stopf('derived_phys')

!$omp parallel do private (lchnk, ncols, ilyr, icol)
   do lchnk = begchunk, endchunk
      ncols=get_ncols_p(lchnk)
      if (pcols > ncols) then
         phys_state(lchnk)%phis(ncols+1:) = 0.0_r8
      end if
   end do
end subroutine d_p_coupling

!=========================================================================================

subroutine p_d_coupling(phys_state, phys_tend, dyn_in, tl_f, tl_qdp)

   use bndry_mod,              only: bndry_exchange
   use edge_mod,               only: edgeVpack, edgeVunpack
   use fvm_mapping,            only: phys2dyn_forcings_fvm
   use test_fvm_mapping,       only: test_mapping_overwrite_tendencies
   use test_fvm_mapping,       only: test_mapping_output_mapped_tendencies
   ! arguments
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend),  intent(inout), dimension(begchunk:endchunk) :: phys_tend
   integer, intent(in)                                              :: tl_qdp, tl_f
   type(dyn_import_t),  intent(inout)                               :: dyn_in
   type(hybrid_t)                                                   :: hybrid

   ! LOCAL VARIABLES
   integer                      :: ic , ncols             ! index
   type(element_t), pointer     :: elem(:)                ! pointer to dyn_in element array
   integer                      :: ie, iep                ! indices over elements
   integer                      :: lchnk, icol, ilyr      ! indices over chunks, columns, layers

   real (kind=r8),  allocatable :: dp_phys(:,:,:)         ! temp array to hold dp on physics grid
   real (kind=r8),  allocatable :: T_tmp(:,:,:)           ! temp array to hold T
   real (kind=r8),  allocatable :: dq_tmp(:,:,:,:)        ! temp array to hold q
   real (kind=r8),  allocatable :: uv_tmp(:,:,:,:)        ! temp array to hold uv
   integer                      :: ioff, m, i, j, k
   integer                      :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)

   integer                      :: tsize                  ! amount of data per grid point passed to physics
   integer                      :: cpter(pcols,0:pver)    ! offsets into chunk buffer for packing data
   integer,         allocatable :: bpter(:,:)             ! offsets into block buffer for unpacking data

   real (kind=r8),  allocatable :: bbuffer(:), cbuffer(:) ! transpose buffers

   real (kind=r8)               :: factor

   integer                      :: num_trac
   integer                      :: nets,nete
   integer                      :: kptr,ii
   !----------------------------------------------------------------------------
   if (iam < par%nprocs) then
      elem => dyn_in%elem
   else
      nullify(elem)
   end if

   allocate(T_tmp(nphys_pts,pver,nelemd))
   allocate(uv_tmp(nphys_pts,2,pver,nelemd))
   allocate(dq_tmp(nphys_pts,pver,pcnst,nelemd))
   allocate(dp_phys(nphys_pts,pver,nelemd))

   T_tmp  = 0.0_r8
   uv_tmp = 0.0_r8
   dq_tmp = 0.0_r8

   if (.not. allocated(q_prev)) then
      call endrun('p_d_coupling: q_prev not allocated')
   end if

   call t_startf('pd_copy')
   if (local_dp_map) then

      !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, ilyr, m, factor)
      do lchnk = begchunk, endchunk
         ncols = get_ncols_p(lchnk)
         call get_gcol_all_p(lchnk, pcols, pgcols)

         call test_mapping_overwrite_tendencies(phys_state(lchnk), phys_tend(lchnk), ncols,&
            lchnk, q_prev(1:ncols,:,:,lchnk))

         do icol = 1, ncols
            call get_gcol_block_d(pgcols(icol), 1, idmb1, idmb2, idmb3)
            ie   = idmb3(1)
            ioff = idmb2(1)

            do ilyr = 1, pver

               ! phys_state%psdry units is Pa
               dp_phys(ioff,ilyr,ie)  = phys_state(lchnk)%pdeldry(icol,ilyr)
               T_tmp(ioff,ilyr,ie)    = phys_tend(lchnk)%dtdt(icol,ilyr)
               uv_tmp(ioff,1,ilyr,ie) = phys_tend(lchnk)%dudt(icol,ilyr)
               uv_tmp(ioff,2,ilyr,ie) = phys_tend(lchnk)%dvdt(icol,ilyr)

               ! convert wet mixing ratios to dry
               ! this code is equivalent to dme_adjust
               factor = phys_state(lchnk)%pdel(icol,ilyr)/phys_state(lchnk)%pdeldry(icol,ilyr)
               do m = 1, pcnst
                  if (cnst_type(m) == 'wet') then
                     phys_state(lchnk)%q(icol,ilyr,m) = factor*phys_state(lchnk)%q(icol,ilyr,m)
                  end if
                  dq_tmp(ioff,ilyr,m,ie) = (phys_state(lchnk)%q(icol,ilyr,m) - &
                                            q_prev(icol,ilyr,m,lchnk))
               end do
            end do
         end do
      end do

   else  ! not local map

      tsize = 4 + pcnst

      allocate(bbuffer(tsize*block_buf_nrecs))
      allocate(cbuffer(tsize*chunk_buf_nrecs))

      !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, cpter, i, icol, ilyr, m, factor)
      do lchnk = begchunk, endchunk
         ncols = get_ncols_p(lchnk)

         call test_mapping_overwrite_tendencies(phys_state(lchnk), phys_tend(lchnk), ncols, lchnk, &
                                                q_prev(1:ncols,:,:,lchnk))

         call chunk_to_block_send_pters(lchnk, pcols, pver+1, tsize, cpter)

         do i = 1, ncols
            cbuffer(cpter(i,0):cpter(i,0)+2+pcnst) = 0.0_r8
         end do

         do icol = 1, ncols
            do ilyr = 1, pver
               cbuffer(cpter(icol,ilyr))   = phys_tend(lchnk)%dtdt(icol,ilyr)
               cbuffer(cpter(icol,ilyr)+1) = phys_tend(lchnk)%dudt(icol,ilyr)
               cbuffer(cpter(icol,ilyr)+2) = phys_tend(lchnk)%dvdt(icol,ilyr)
               cbuffer(cpter(icol,ilyr)+3) = phys_state(lchnk)%pdeldry(icol,ilyr)

               ! this code is equivalent to dme_adjust
               factor = phys_state(lchnk)%pdel(icol,ilyr)/phys_state(lchnk)%pdeldry(icol,ilyr)

               do m = 1, pcnst
                  if (cnst_type(m) == 'wet') then
                     phys_state(lchnk)%q(icol,ilyr,m) = factor*phys_state(lchnk)%q(icol,ilyr,m)
                  end if
                  cbuffer(cpter(icol,ilyr)+3+m) = (phys_state(lchnk)%q(icol,ilyr,m) - &
                                                   q_prev(icol,ilyr,m,lchnk))
               end do
            end do
          end do
      end do

      call t_barrierf('sync_chk_to_blk', mpicom)
      call t_startf ('chunk_to_block')
      call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
      call t_stopf  ('chunk_to_block')

      if (iam < par%nprocs) then

         if (fv_nphys > 0) then
            allocate(bpter(fv_nphys*fv_nphys,0:pver))
         else
            allocate(bpter(npsq,0:pver))
         end if

         !$omp parallel do num_threads(horz_num_threads) private (ie, bpter, icol, ilyr, m, ncols)
         do ie = 1, nelemd

            if (fv_nphys > 0) then
               call chunk_to_block_recv_pters(elem(ie)%GlobalID, fv_nphys*fv_nphys, &
                                              pver+1, tsize, bpter)
               ncols = fv_nphys*fv_nphys
            else
               call chunk_to_block_recv_pters(elem(ie)%GlobalID, npsq, &
                                              pver+1, tsize, bpter)
               ncols = elem(ie)%idxP%NumUniquePts
            end if

            do icol = 1, ncols
               do ilyr = 1, pver
                  T_tmp   (icol,ilyr,ie)   = bbuffer(bpter(icol,ilyr))
                  uv_tmp  (icol,1,ilyr,ie) = bbuffer(bpter(icol,ilyr)+1)
                  uv_tmp  (icol,2,ilyr,ie) = bbuffer(bpter(icol,ilyr)+2)
                  dp_phys (icol,ilyr,ie)   = bbuffer(bpter(icol,ilyr)+3)

                  do m = 1, pcnst
                     dq_tmp(icol,ilyr,m,ie) = bbuffer(bpter(icol,ilyr)+3+m)
                  end do
               end do
            end do
         end do
         deallocate(bpter)

      end if

      deallocate( bbuffer )
      deallocate( cbuffer )

   end if
   call t_stopf('pd_copy')

   if (iam < par%nprocs) then

      if (fv_nphys > 0) then

         ! put forcings into fvm structure
         num_trac = max(qsize,ntrac)
         do ie = 1, nelemd
            do j = 1, fv_nphys
               do i = 1, fv_nphys
                  ii = i + (j-1)*fv_nphys
                  dyn_in%fvm(ie)%ft(i,j,1:pver)     = T_tmp(ii,1:pver,ie)
                  dyn_in%fvm(ie)%fm(i,j,1:2,1:pver) = uv_tmp(ii,1:2,1:pver,ie)
                  dyn_in%fvm(ie)%fc_phys(i,j,1:pver,1:num_trac) = dq_tmp(ii,1:pver,1:num_trac,ie)
                  dyn_in%fvm(ie)%dp_phys(i,j,1:pver) = dp_phys(ii,1:pver,ie)
               end do
            end do
         end do

         !JMD $OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete,n)
         !JMD        hybrid = config_thread_region(par,'horizontal')
         hybrid = config_thread_region(par,'serial')
         call get_loop_ranges(hybrid,ibeg=nets,iend=nete)

         ! high-order mapping of ft and fm (and fq if no cslam) using fvm technology
         call t_startf('phys2dyn')
         call phys2dyn_forcings_fvm(elem, dyn_in%fvm, hybrid,nets,nete,ntrac==0, tl_f, tl_qdp)
         call t_stopf('phys2dyn')
      else

         call t_startf('putUniquePoints')

         !$omp parallel do num_threads(horz_num_threads) private(ie,ncols)
         do ie = 1, nelemd
            ncols = elem(ie)%idxP%NumUniquePts
            call putUniquePoints(elem(ie)%idxP, nlev, T_tmp(1:ncols,:,ie),       &
               elem(ie)%derived%fT(:,:,:))
            call putUniquePoints(elem(ie)%idxV, 2, nlev, uv_tmp(1:ncols,:,:,ie), &
               elem(ie)%derived%fM(:,:,:,:))
            call putUniquePoints(elem(ie)%idxV, nlev,pcnst, dq_tmp(1:ncols,:,:,ie), &
               elem(ie)%derived%fQ(:,:,:,:))
         end do
         call t_stopf('putUniquePoints')
      end if
   end if

   deallocate(T_tmp)
   deallocate(uv_tmp)
   deallocate(dq_tmp)

   ! Boundary exchange for physics forcing terms.
   ! For physics on GLL grid, for points with duplicate degrees of freedom,
   ! putuniquepoints() set one of the element values and set the others to zero,
   ! so do a simple sum (boundary exchange with no weights).
   ! For physics grid, we interpolated into all points, so do weighted average.

   call t_startf('bndry_exchange')

   do ie = 1, nelemd
      if (fv_nphys > 0) then
         do k = 1, nlev
            dyn_in%elem(ie)%derived%FM(:,:,1,k) =                          &
                 dyn_in%elem(ie)%derived%FM(:,:,1,k) *                     &
                 dyn_in%elem(ie)%spheremp(:,:)
            dyn_in%elem(ie)%derived%FM(:,:,2,k) =                          &
                 dyn_in%elem(ie)%derived%FM(:,:,2,k) *                     &
                 dyn_in%elem(ie)%spheremp(:,:)
            dyn_in%elem(ie)%derived%FT(:,:,k) =                            &
                 dyn_in%elem(ie)%derived%FT(:,:,k) *                       &
                 dyn_in%elem(ie)%spheremp(:,:)
            do m = 1, qsize
               dyn_in%elem(ie)%derived%FQ(:,:,k,m) =                       &
                    dyn_in%elem(ie)%derived%FQ(:,:,k,m) *                  &
                    dyn_in%elem(ie)%spheremp(:,:)
            end do
         end do
      end if
      kptr = 0
      call edgeVpack(edgebuf, dyn_in%elem(ie)%derived%FM(:,:,:,:), 2*nlev, kptr, ie)
      kptr = kptr + 2*nlev
      call edgeVpack(edgebuf, dyn_in%elem(ie)%derived%FT(:,:,:), nlev, kptr, ie)
      kptr = kptr + nlev
      call edgeVpack(edgebuf, dyn_in%elem(ie)%derived%FQ(:,:,:,:), nlev*qsize, kptr, ie)
   end do

   if (iam < par%nprocs) then
     call bndry_exchange(par, edgebuf, location='p_d_coupling')
   end if

   do ie = 1, nelemd
      kptr = 0
      call edgeVunpack(edgebuf, dyn_in%elem(ie)%derived%FM(:,:,:,:), 2*nlev, kptr, ie)
      kptr = kptr + 2*nlev
      call edgeVunpack(edgebuf, dyn_in%elem(ie)%derived%FT(:,:,:), nlev, kptr, ie)
      kptr = kptr + nlev
      call edgeVunpack(edgebuf, dyn_in%elem(ie)%derived%FQ(:,:,:,:), nlev*qsize, kptr, ie)
      if (fv_nphys > 0) then
         do k = 1, nlev
            dyn_in%elem(ie)%derived%FM(:,:,1,k) =                             &
                 dyn_in%elem(ie)%derived%FM(:,:,1,k) *                        &
                 dyn_in%elem(ie)%rspheremp(:,:)
            dyn_in%elem(ie)%derived%FM(:,:,2,k) =                             &
                 dyn_in%elem(ie)%derived%FM(:,:,2,k) *                        &
                 dyn_in%elem(ie)%rspheremp(:,:)
            dyn_in%elem(ie)%derived%FT(:,:,k) =                               &
                 dyn_in%elem(ie)%derived%FT(:,:,k) *                          &
                 dyn_in%elem(ie)%rspheremp(:,:)
            do m = 1, qsize
               dyn_in%elem(ie)%derived%FQ(:,:,k,m) =                          &
                    dyn_in%elem(ie)%derived%FQ(:,:,k,m) *                     &
                    dyn_in%elem(ie)%rspheremp(:,:)
            end do
         end do
      end if
   end do
   call t_stopf('bndry_exchange')

   if (iam < par%nprocs .and. fv_nphys > 0) then
      call test_mapping_output_mapped_tendencies(dyn_in%fvm(1:nelemd), elem(1:nelemd), &
                                                 1, nelemd, tl_f, tl_qdp)
   end if
end subroutine p_d_coupling

!=========================================================================================

subroutine derived_phys_dry(phys_state, phys_tend, pbuf2d)

   use constituents,  only: qmin
   use physconst,     only: cpair, gravit, rair, zvir, cappa, rairv,rh2o,rair
   use geopotential,  only: geopotential_t
   use physics_types, only: set_state_pdry, set_wet_to_dry
   use check_energy,  only: check_energy_timestep_init
   use hycoef,        only: hyam, hybm, hyai, hybi, ps0
   use shr_vmath_mod, only: shr_vmath_log
   !use gmean_mod,     only: gmean
   use phys_gmean,    only: gmean
   use qneg_module,   only: qneg3_new

   ! arguments
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
   type(physics_buffer_desc),      pointer     :: pbuf2d(:,:)

   ! local variables
   integer :: lchnk
   real(r8) :: qbot                 ! bottom level q before change
   real(r8) :: qbotm1               ! bottom-1 level q before change
   real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
   real(r8) :: qmavl                ! available q at level pver-1

   real(r8) :: ke(pcols,begchunk:endchunk)
   real(r8) :: se(pcols,begchunk:endchunk)
   real(r8) :: ke_glob(1),se_glob(1)
   real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer
   real(r8) :: factor_array(pcols,nlev)

   integer :: m, i, k, ncol

   type(physics_buffer_desc), pointer :: pbuf_chnk(:)
   !----------------------------------------------------------------------------

   ! Evaluate derived quantities

   !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncol, k, i, m , zvirv, pbuf_chnk, factor_array)
   do lchnk = begchunk,endchunk

      ncol = get_ncols_p(lchnk)

      ! dry pressure variables

      do i = 1, ncol
         phys_state(lchnk)%psdry(i) = hyai(1)*ps0 + sum(phys_state(lchnk)%pdel(i,:))
      end do
      do i = 1, ncol
         phys_state(lchnk)%pintdry(i,1) = hyai(1)*ps0
      end do
      call shr_vmath_log(phys_state(lchnk)%pintdry(1:ncol,1), &
                         phys_state(lchnk)%lnpintdry(1:ncol,1),ncol)
      do k = 1, nlev
         do i = 1, ncol
            phys_state(lchnk)%pintdry(i,k+1) = phys_state(lchnk)%pintdry(i,k) + &
                                               phys_state(lchnk)%pdel(i,k)
         end do
         call shr_vmath_log(phys_state(lchnk)%pintdry(1:ncol,k+1),&
                            phys_state(lchnk)%lnpintdry(1:ncol,k+1),ncol)
      end do

      do k=1,nlev
         do i=1,ncol
            phys_state(lchnk)%pdeldry (i,k) = phys_state(lchnk)%pdel(i,k)
            phys_state(lchnk)%rpdeldry(i,k) = 1._r8/phys_state(lchnk)%pdeldry(i,k)
            phys_state(lchnk)%pmiddry (i,k) = 0.5D0*(phys_state(lchnk)%pintdry(i,k+1) + &
                                                     phys_state(lchnk)%pintdry(i,k))
         end do
         call shr_vmath_log(phys_state(lchnk)%pmiddry(1:ncol,k), &
                            phys_state(lchnk)%lnpmiddry(1:ncol,k),ncol)
      end do

      ! wet pressure variables (should be removed from physics!)

      do k=1,nlev
         do i=1,ncol
            ! to be consistent with total energy formula in physic's check_energy module only
            ! include water vapor in in moist dp
            factor_array(i,k) = 1+phys_state(lchnk)%q(i,k,1)
         end do
      end do

      do k=1,nlev
         do i=1,ncol
            phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pdeldry(i,k)*factor_array(i,k)
         end do
      end do

      ! initialize vertical loop - model top pressure

      do i=1,ncol
         phys_state(lchnk)%ps(i)     = phys_state(lchnk)%pintdry(i,1)
         phys_state(lchnk)%pint(i,1) = phys_state(lchnk)%pintdry(i,1)
      end do
      do k = 1, nlev
         do i=1,ncol
            phys_state(lchnk)%pint(i,k+1) =  phys_state(lchnk)%pint(i,k)+phys_state(lchnk)%pdel(i,k)
            phys_state(lchnk)%pmid(i,k)   = (phys_state(lchnk)%pint(i,k+1)+phys_state(lchnk)%pint(i,k))/2._r8
            phys_state(lchnk)%ps  (i)     =  phys_state(lchnk)%ps(i) + phys_state(lchnk)%pdel(i,k)
         end do
         call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k),phys_state(lchnk)%lnpint(1:ncol,k),ncol)
         call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k),phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
      end do
      call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp),phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)

      do k = 1, nlev
         do i = 1, ncol
            phys_state(lchnk)%rpdel(i,k)  = 1._r8/phys_state(lchnk)%pdel(i,k)
            phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                                            / phys_state(lchnk)%pmid(i,k))**cappa
         end do
      end do

      ! all tracers (including moisture) are in dry mixing ratio units
      ! physics expect water variables moist
      factor_array(1:ncol,1:nlev) = 1/factor_array(1:ncol,1:nlev)

      do m = 1,pcnst
         if (cnst_type(m) == 'wet') then
            do k = 1, nlev
               do i = 1, ncol
                  phys_state(lchnk)%q(i,k,m) = factor_array(i,k)*phys_state(lchnk)%q(i,k,m)
               end do
            end do
         end if
      end do


      ! fill zvirv 2D variables to be compatible with geopotential_t interface
      zvirv(:,:) = zvir

      ! Compute initial geopotential heights - based on full pressure
      call geopotential_t (phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
         phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
         phys_state(lchnk)%t     , phys_state(lchnk)%q(:,:,1), rairv(:,:,lchnk),  gravit,  zvirv       , &
         phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol                )

      ! Compute initial dry static energy, include surface geopotential
      do k = 1, pver
         do i = 1, ncol
#if FIX_TOTE
            ! general formula:  E = CV_air T + phis + gravit*zi )
            ! hydrostatic case: integrate zi term by parts, use CP=CV+R to get:
            ! E = CP_air T + phis   (Holton Section 8.3)
            ! to use this, update geopotential.F90, and other not-yet-found physics routines:
            ! (check boundary layer code, others which have gravit and zi() or zm()
            phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
               + phys_state(lchnk)%phis(i)
#else
            phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
               + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
#endif
         end do
      end do

      ! Ensure tracers are all positive
      call qneg3_new('D_P_COUPLING',lchnk  ,ncol    ,pcols   ,pver    , &
           1, pcnst, qmin  ,phys_state(lchnk)%q)

      ! Compute energy and water integrals of input state
      pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
      call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)


   end do  ! lchnk

end subroutine derived_phys_dry

end module dp_coupling
