!---------------------------------------------------------------------------------------------------
! dynamics - physics coupling module
!---------------------------------------------------------------------------------------------------
module dp_coupling
  use constituents,   only: pcnst, cnst_name
  use cam_history,    only: outfld, write_inithist, hist_fld_active
  use dimensions_mod, only: np, npsq, nelemd, nlev
  use dof_mod,        only: UniquePoints, PutUniquePoints
  use dyn_comp,       only: dyn_export_t, dyn_import_t, TimeLevel
  use dyn_grid,       only: get_gcol_block_d, fv_nphys
  use pmgrid,         only: plev
  use element_mod,    only: element_t
  use kinds,          only: real_kind, int_kind
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use physics_types,  only: physics_state, physics_tend 
  use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
  use cam_logfile,    only: iulog
  use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
  use spmd_utils,     only: mpicom, iam
  use perf_mod,       only: t_startf, t_stopf, t_barrierf
  use parallel_mod,   only: par
  use iop_data_mod,   only: single_column
  use element_ops,    only: get_temperature
  use phys_grid,      only: get_ncols_p, get_gcol_all_p, &
                            transpose_block_to_chunk, transpose_chunk_to_block,   &
                            chunk_to_block_send_pters, chunk_to_block_recv_pters, &
                            block_to_chunk_recv_pters, block_to_chunk_send_pters
  private
  public :: d_p_coupling, p_d_coupling

CONTAINS
  !=================================================================================================
  !=================================================================================================
  subroutine d_p_coupling(phys_state, phys_tend,  pbuf2d, dyn_out)
    use physics_buffer,          only: physics_buffer_desc, pbuf_get_chunk, pbuf_get_field
    use shr_vmath_mod,           only: shr_vmath_exp
    use time_manager,            only: is_first_step
    use cam_abortutils,          only: endrun
    use gravity_waves_sources,   only: gws_src_fnct
    use dyn_comp,                only: frontgf_idx, frontga_idx, hvcoord
    use phys_control,            only: use_gw_front
    use dyn_comp,                only: dom_mt
    use gllfvremap_mod,          only: gfr_dyn_to_fv_phys

    implicit none
    !---------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    type(dyn_export_t), intent(inout)  :: dyn_out         ! dynamics export 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)     ! physics buffer
    ! OUTPUT PARAMETERS:
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend 
    ! LOCAL VARIABLES
    real(kind=real_kind), dimension(npsq,nelemd)            :: ps_tmp ! temp array to hold ps
    real(kind=real_kind), dimension(npsq,nelemd)            :: zs_tmp ! temp array to hold phis  
    real(kind=real_kind), dimension(npsq,pver,nelemd)       :: T_tmp  ! temp array to hold T
    real(kind=real_kind), dimension(npsq,2,pver,nelemd)     :: uv_tmp ! temp array to hold u and v
    real(kind=real_kind), dimension(npsq,pver,pcnst,nelemd) :: q_tmp  ! temp to hold advected constituents
    real(kind=real_kind), dimension(npsq,pver,nelemd)       :: om_tmp ! temp array to hold omega
    type(element_t),          pointer :: elem(:)          ! pointer to dyn_out element array
    type(physics_buffer_desc),pointer :: pbuf_chnk(:)     ! temporary pbuf pointer
    integer(kind=int_kind)   :: ie                        ! indices over elements
    integer(kind=int_kind)   :: lchnk, icol, ilyr         ! indices over chunks, columns, layers
    integer                  :: tl_f                      ! time level
    integer                  :: ncols, ierr, ncol_d       ! 
    integer                  :: i,j,m                     ! loop iterators
    integer(kind=int_kind)   :: ioff                      ! column index within block
    integer                  :: pgcols(pcols)             ! global column indices for chunk
    integer                  :: idmb1(1)                  ! block id
    integer                  :: idmb2(1)                  ! column index within block
    integer                  :: idmb3(1)                  ! local block id
    integer                  :: tsize                     ! # vars per grid point passed to physics
    integer                  :: bpter(npsq,0:pver)        ! offsets into block buffer for packing 
    integer                  :: cpter(pcols,0:pver)       ! offsets into chunk buffer for unpacking 
    integer                  :: nphys, nphys_sq           ! physics grid parameters
    real (kind=real_kind)    :: temperature(np,np,nlev)   ! Temperature from dynamics
    ! Frontogenesis
    real (kind=real_kind), allocatable :: frontgf(:,:,:)  ! frontogenesis function
    real (kind=real_kind), allocatable :: frontga(:,:,:)  ! frontogenesis angle 
    real (kind=r8),            pointer :: pbuf_frontgf(:,:)
    real (kind=r8),            pointer :: pbuf_frontga(:,:)
    ! Transpose buffers
    real (kind=real_kind), allocatable, dimension(:) :: bbuffer 
    real (kind=real_kind), allocatable, dimension(:) :: cbuffer
    !---------------------------------------------------------------------------

    nullify(pbuf_chnk)
    nullify(pbuf_frontgf)
    nullify(pbuf_frontga)

    if (fv_nphys > 0) then
      nphys = fv_nphys
    else
      nphys = np
    end if
    nphys_sq = nphys*nphys

    if (use_gw_front) then
       allocate(frontgf(nphys_sq,pver,nelemd), stat=ierr)
       if (ierr /= 0) call endrun("dp_coupling: Allocate of frontgf failed.")
       allocate(frontga(nphys_sq,pver,nelemd), stat=ierr)
       if (ierr /= 0) call endrun("dp_coupling: Allocate of frontga failed.")
    end if

    if( par%dynproc) then

      elem => dyn_out%elem
      tl_f = TimeLevel%n0  ! time split physics (with forward-in-time RK)

      if (use_gw_front) call gws_src_fnct(elem, tl_f, nphys, frontgf, frontga)

      if (fv_nphys > 0) then
        !-----------------------------------------------------------------------
        ! Map dynamics state to FV physics grid
        !-----------------------------------------------------------------------
        call t_startf('dyn_to_fv_phys')
        call gfr_dyn_to_fv_phys(par, dom_mt, tl_f, hvcoord, elem, ps_tmp, zs_tmp, &
             T_tmp, uv_tmp, om_tmp, q_tmp)
        call t_stopf('dyn_to_fv_phys')

        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
      else ! fv_nphys > 0
        !-----------------------------------------------------------------------
        ! Physics on GLL grid: collect unique points before copying
        !-----------------------------------------------------------------------
        call t_startf('UniquePoints')
        do ie = 1,nelemd
          ncols = elem(ie)%idxP%NumUniquePts
          call get_temperature(elem(ie),temperature,hvcoord,tl_f)
          call UniquePoints(elem(ie)%idxP,       elem(ie)%state%ps_v(:,:,tl_f), ps_tmp(1:ncols,ie))
          call UniquePoints(elem(ie)%idxP,       elem(ie)%state%phis,           zs_tmp(1:ncols,ie))
          call UniquePoints(elem(ie)%idxP,  nlev, temperature,                  T_tmp(1:ncols,:,ie))
          call UniquePoints(elem(ie)%idxP,  nlev,elem(ie)%derived%omega_p,      om_tmp(1:ncols,:,ie))
          call UniquePoints(elem(ie)%idxP,2,nlev,elem(ie)%state%V(:,:,:,:,tl_f),uv_tmp(1:ncols,:,:,ie))
          call UniquePoints(elem(ie)%idxP,nlev,pcnst,elem(ie)%state%Q(:,:,:,:), q_tmp(1:ncols,:,:,ie))
        end do
        call t_stopf('UniquePoints')
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
      end if ! fv_nphys > 0

    else ! par%dynproc

      ps_tmp(:,:)      = 0._r8
      T_tmp(:,:,:)     = 0._r8
      uv_tmp(:,:,:,:)  = 0._r8
      om_tmp(:,:,:)    = 0._r8
      zs_tmp(:,:)      = 0._r8
      q_tmp(:,:,:,:)   = 0._r8
      if (use_gw_front) then
        frontgf(:,:,:) = 0._r8
        frontga(:,:,:) = 0._r8
      end if

    end if ! par%dynproc

    call t_startf('dpcopy')
    if (local_dp_map) then

      !$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, ilyr, m, pbuf_chnk, pbuf_frontgf, pbuf_frontga)
      do lchnk = begchunk,endchunk
        ncols=get_ncols_p(lchnk)
        call get_gcol_all_p(lchnk,pcols,pgcols)

        pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

        if (use_gw_front) then
          call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
          call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
        end if

        do icol = 1,ncols
          call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
          ie = idmb3(1)
          ioff = idmb2(1)
          phys_state(lchnk)%ps(icol)   = ps_tmp(ioff,ie)
          phys_state(lchnk)%phis(icol) = zs_tmp(ioff,ie)
          do ilyr = 1,pver
            phys_state(lchnk)%t(icol,ilyr)     = T_tmp(ioff,ilyr,ie)	   
            phys_state(lchnk)%u(icol,ilyr)     = uv_tmp(ioff,1,ilyr,ie)
            phys_state(lchnk)%v(icol,ilyr)     = uv_tmp(ioff,2,ilyr,ie)
            phys_state(lchnk)%omega(icol,ilyr) = om_tmp(ioff,ilyr,ie)
            if (use_gw_front) then
              pbuf_frontgf(icol,ilyr) = frontgf(ioff,ilyr,ie)
              pbuf_frontga(icol,ilyr) = frontga(ioff,ilyr,ie)
            end if
          end do ! ilyr

          do m = 1,pcnst
            do ilyr = 1,pver
              phys_state(lchnk)%q(icol,ilyr,m) = q_tmp(ioff,ilyr,m,ie)
            end do ! ilyr
          end do ! m
        end do ! icol
      end do ! lchnk

    else  ! .not. local_dp_map

      tsize = 4 + pcnst
      if (use_gw_front) tsize = tsize + 2

      allocate(bbuffer(tsize*block_buf_nrecs))
      allocate(cbuffer(tsize*chunk_buf_nrecs))

      if (par%dynproc) then

        !$omp parallel do private (ie, bpter, icol, ilyr, m, ncols)
        do ie = 1,nelemd
          call block_to_chunk_send_pters(elem(ie)%GlobalID,nphys_sq,pver+1,tsize,bpter(1:nphys_sq,:))
          if (fv_nphys > 0) then
            ncols = nphys_sq
          else
            ncols = elem(ie)%idxP%NumUniquePts
          end if
          do icol = 1,ncols
            bbuffer(bpter(icol,0)+2:bpter(icol,0)+tsize-1) = 0.0_r8
            bbuffer(bpter(icol,0))   = ps_tmp(icol,ie)
            bbuffer(bpter(icol,0)+1) = zs_tmp(icol,ie)
            do ilyr = 1,pver
              bbuffer(bpter(icol,ilyr))   = T_tmp(icol,ilyr,ie)
              bbuffer(bpter(icol,ilyr)+1) = uv_tmp(icol,1,ilyr,ie)
              bbuffer(bpter(icol,ilyr)+2) = uv_tmp(icol,2,ilyr,ie)
              bbuffer(bpter(icol,ilyr)+3) = om_tmp(icol,ilyr,ie)
              if (use_gw_front) then
                bbuffer(bpter(icol,ilyr)+4) = frontgf(icol,ilyr,ie)
                bbuffer(bpter(icol,ilyr)+5) = frontga(icol,ilyr,ie)
              end if
              do m = 1,pcnst
                bbuffer(bpter(icol,ilyr)+tsize-pcnst-1+m) = q_tmp(icol,ilyr,m,ie)
              end do
            end do ! ilyr
          end do ! icol
        end do ! ie

      else
        bbuffer(:) = 0._r8
      end if ! par%dynproc

      call t_barrierf ('sync_blk_to_chk', mpicom)
      call t_startf ('block_to_chunk')
      call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
      call t_stopf  ('block_to_chunk')

      !$omp parallel do private (lchnk, ncols, cpter, icol, ilyr, m, pbuf_chnk, pbuf_frontgf, pbuf_frontga)
      do lchnk = begchunk,endchunk
        ncols = phys_state(lchnk)%ncol
        pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
        if (use_gw_front) then
           call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
           call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
        end if
        call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)
        do icol = 1,ncols
          phys_state(lchnk)%ps  (icol) = cbuffer(cpter(icol,0))
          phys_state(lchnk)%phis(icol) = cbuffer(cpter(icol,0)+1)
          do ilyr = 1,pver
            phys_state(lchnk)%t    (icol,ilyr) = cbuffer(cpter(icol,ilyr))
            phys_state(lchnk)%u    (icol,ilyr) = cbuffer(cpter(icol,ilyr)+1)
            phys_state(lchnk)%v    (icol,ilyr) = cbuffer(cpter(icol,ilyr)+2)
            phys_state(lchnk)%omega(icol,ilyr) = cbuffer(cpter(icol,ilyr)+3)
             if (use_gw_front) then
                pbuf_frontgf(icol,ilyr) = cbuffer(cpter(icol,ilyr)+4)
                pbuf_frontga(icol,ilyr) = cbuffer(cpter(icol,ilyr)+5)
             end if
             do m = 1,pcnst
                phys_state(lchnk)%q(icol,ilyr,m) = cbuffer(cpter(icol,ilyr)+tsize-pcnst-1+m)
             end do ! m
          end do ! ilyr
        end do ! icol
      end do ! lchnk

      deallocate( bbuffer )
      deallocate( cbuffer )

    end if ! local_dp_map
    call t_stopf('dpcopy')

    call t_startf('derived_phys')
    call derived_phys(phys_state,phys_tend,pbuf2d)
    call t_stopf('derived_phys')

!for theta there is no need to multiply omega_p by p
#ifndef MODEL_THETA_L
    !$omp parallel do private (lchnk, ncols, ilyr, icol)
    do lchnk = begchunk,endchunk
      ncols = get_ncols_p(lchnk)
      do ilyr = 1,pver
        do icol = 1,ncols
          if (.not.single_column) then
            phys_state(lchnk)%omega(icol,ilyr) = phys_state(lchnk)%omega(icol,ilyr) &
                                                *phys_state(lchnk)%pmid(icol,ilyr)
          end if
        end do ! icol
      end do ! ilyr
    end do ! lchnk
#endif

    if ( write_inithist() ) then
      if (fv_nphys > 0) then

        ncol_d = np*np
        do ie = 1,nelemd
          ncols = elem(ie)%idxP%NumUniquePts
          call outfld('PS&IC',elem(ie)%state%ps_v(:,:,tl_f),  ncol_d,ie)
          call outfld('U&IC', elem(ie)%state%V(:,:,1,:,tl_f), ncol_d,ie)
          call outfld('V&IC', elem(ie)%state%V(:,:,2,:,tl_f), ncol_d,ie)
          call get_temperature(elem(ie),temperature,hvcoord,tl_f)
          call outfld('T&IC',temperature,ncol_d,ie)
          do m = 1,pcnst
            call outfld(trim(cnst_name(m))//'&IC',elem(ie)%state%Q(:,:,:,m), ncol_d,ie)
          end do ! m
        end do ! ie
        
      else

        do lchnk = begchunk,endchunk
          call outfld('T&IC', phys_state(lchnk)%t, pcols,lchnk)
          call outfld('U&IC', phys_state(lchnk)%u, pcols,lchnk)
          call outfld('V&IC', phys_state(lchnk)%v, pcols,lchnk)
          call outfld('PS&IC',phys_state(lchnk)%ps,pcols,lchnk)
          do m = 1,pcnst
            call outfld(trim(cnst_name(m))//'&IC',phys_state(lchnk)%q(1,1,m), pcols,lchnk)
          end do ! m
        end do ! lchnk

      end if ! fv_nphys > 0
    end if ! write_inithist
   
  end subroutine d_p_coupling
  !=================================================================================================
  !=================================================================================================
  subroutine p_d_coupling(phys_state, phys_tend,  dyn_in)
    use shr_vmath_mod,           only: shr_vmath_log
    use cam_control_mod,         only: adiabatic
    use control_mod,             only: ftype
    use dyn_comp,                only: dom_mt, hvcoord
    use gllfvremap_mod,          only: gfr_fv_phys_to_dyn
    use time_manager,            only: get_step_size
    implicit none
    ! INPUT PARAMETERS:
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend),  intent(inout), dimension(begchunk:endchunk) :: phys_tend
    ! OUTPUT PARAMETERS:
    type(dyn_import_t),  intent(inout)   :: dyn_in
    ! LOCAL VARIABLES
    integer :: ic , ncols                                  ! index
    type(element_t), pointer :: elem(:)                    ! pointer to dyn_in element array
    integer(kind=int_kind)   :: ie, iep                    ! indices over elements
    integer(kind=int_kind)   :: lchnk, icol, ilyr          ! indices over chunks, columns, layers
    real (kind=real_kind), dimension(npsq,pver,nelemd)       :: T_tmp  ! temp array to hold T
    real (kind=real_kind), dimension(npsq,2,pver,nelemd)     :: uv_tmp ! temp array to hold u and v
    real (kind=real_kind), dimension(npsq,pver,pcnst,nelemd) :: q_tmp  ! temp to hold advected constituents
    real (kind=real_kind)    :: dtime
    integer(kind=int_kind)   :: m, i, j, k                 ! loop iterators
    integer(kind=int_kind)   :: gi(2), gj(2)               ! index list used to simplify pg2 case
    integer(kind=int_kind)   :: di, dj
    integer(kind=int_kind)   :: pgcols(pcols)
    integer(kind=int_kind)   :: ioff
    integer(kind=int_kind)   :: idmb1(1)
    integer(kind=int_kind)   :: idmb2(1)
    integer(kind=int_kind)   :: idmb3(1)
    integer                  :: nphys, nphys_sq
    integer                  :: tsize                 ! # vars per grid point passed to physics
    integer                  :: cpter(pcols,0:pver)   ! offsets into chunk buffer for packing 
    integer                  :: bpter(npsq,0:pver)    ! offsets into block buffer for unpacking 
    ! Transpose buffers
    real (kind=real_kind), allocatable, dimension(:) :: bbuffer, cbuffer 
    !---------------------------------------------------------------------------
    if (par%dynproc) then
      elem => dyn_in%elem
    else
      nullify(elem)
    end if

    if (fv_nphys > 0) then
      nphys = fv_nphys
    else
      nphys = np
    end if
    nphys_sq = nphys*nphys

    T_tmp  = 0.0_r8
    uv_tmp = 0.0_r8
    q_tmp  = 0.0_r8

    if(adiabatic) return

    call t_startf('pd_copy')
    if(local_dp_map) then

      !$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, ilyr, m)
      do lchnk = begchunk,endchunk
        ncols = get_ncols_p(lchnk)
        call get_gcol_all_p(lchnk,pcols,pgcols)
        do icol = 1,ncols
          call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
          ie = idmb3(1)
          ioff = idmb2(1)
          do ilyr = 1,pver
            T_tmp(ioff,ilyr,ie)    = phys_tend(lchnk)%dtdt(icol,ilyr)
            uv_tmp(ioff,1,ilyr,ie) = phys_tend(lchnk)%dudt(icol,ilyr)
            uv_tmp(ioff,2,ilyr,ie) = phys_tend(lchnk)%dvdt(icol,ilyr)
            do m = 1,pcnst
              q_tmp(ioff,ilyr,m,ie) = phys_state(lchnk)%q(icol,ilyr,m)
            end do
          end do ! ilyr
      	end do ! icol
      end do ! lchnk

    else ! local_dp_map

      tsize = 3 + pcnst

      allocate( bbuffer(tsize*block_buf_nrecs) )
      allocate( cbuffer(tsize*chunk_buf_nrecs) )

      !$omp parallel do private (lchnk, ncols, cpter, i, icol, ilyr, m)
      do lchnk = begchunk,endchunk
        ncols = get_ncols_p(lchnk)
        call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)
        do i = 1,ncols
           cbuffer(cpter(i,0):cpter(i,0)+2+pcnst) = 0.0_r8
        end do
        do icol = 1,ncols
          do ilyr = 1,pver
            cbuffer(cpter(icol,ilyr))   = phys_tend(lchnk)%dtdt(icol,ilyr)
            cbuffer(cpter(icol,ilyr)+1) = phys_tend(lchnk)%dudt(icol,ilyr)
            cbuffer(cpter(icol,ilyr)+2) = phys_tend(lchnk)%dvdt(icol,ilyr)
            do m=1,pcnst
              cbuffer(cpter(icol,ilyr)+2+m) = phys_state(lchnk)%q(icol,ilyr,m)
            end do
          end do ! ilyr
        end do ! icol
      end do ! lchnk

      call t_barrierf('sync_chk_to_blk', mpicom)
      call t_startf ('chunk_to_block')
      call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
      call t_stopf  ('chunk_to_block')

      if (par%dynproc) then
        !$omp parallel do private (ie, bpter, icol, ilyr, m, ncols)
        do ie = 1,nelemd
          if (fv_nphys > 0) then
            ncols = nphys_sq
          else
            ncols = elem(ie)%idxP%NumUniquePts
          end if
          call chunk_to_block_recv_pters(elem(ie)%GlobalID,nphys_sq,pver+1,tsize,bpter(1:nphys_sq,:))
          do icol = 1,ncols
            do ilyr = 1,pver
              T_tmp  (icol,ilyr,ie)   = bbuffer(bpter(icol,ilyr))
              uv_tmp (icol,1,ilyr,ie) = bbuffer(bpter(icol,ilyr)+1)
              uv_tmp (icol,2,ilyr,ie) = bbuffer(bpter(icol,ilyr)+2)
              do m = 1,pcnst
                q_tmp(icol,ilyr,m,ie) = bbuffer(bpter(icol,ilyr)+2+m)
              end do
            end do ! ilyr
          end do ! icol
        end do ! ie
      end if ! par%dynproc

      deallocate( bbuffer )
      deallocate( cbuffer )
       
    end if ! local_dp_map
    call t_stopf('pd_copy')

    if (par%dynproc) then
      if (fv_nphys > 0) then
        call t_startf('fv_phys_to_dyn')
        ! Map FV physics state to dynamics grid
        dtime = get_step_size()
        call gfr_fv_phys_to_dyn(par, dom_mt, TimeLevel%n0, hvcoord, elem, T_tmp, &
             uv_tmp, q_tmp)
        call t_stopf('fv_phys_to_dyn')

      else ! physics is on GLL nodes

        call t_startf('putUniquePoints')
        do ie = 1,nelemd
          ncols = elem(ie)%idxP%NumUniquePts
          call putUniquePoints(elem(ie)%idxP,    nlev,       T_tmp(1:ncols,:,ie),   &
                               elem(ie)%derived%fT(:,:,:))
          call putUniquePoints(elem(ie)%idxP, 2, nlev,       uv_tmp(1:ncols,:,:,ie), &
                               elem(ie)%derived%fM(:,:,:,:))
          call putUniquePoints(elem(ie)%idxP,    nlev,pcnst, q_tmp(1:ncols,:,:,ie),   &
                               elem(ie)%derived%fQ(:,:,:,:))
        end do ! ie
        call t_stopf('putUniquePoints')

      end if ! fv_nphys > 0
    end if ! par%dynproc

  end subroutine p_d_coupling
  !=================================================================================================
  !=================================================================================================
  subroutine derived_phys(phys_state, phys_tend, pbuf2d)
    use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk
    use constituents,   only: qmin
    use physconst,      only: cpair, gravit, rair, zvir, cappa, rairv
    use spmd_utils,     only: masterproc
    use ppgrid,         only: pver
    use geopotential,   only: geopotential_t
    use physics_types,  only: set_state_pdry, set_wet_to_dry
    use check_energy,   only: check_energy_timestep_init
    use hycoef,         only: hyam, hybm, hyai, hybi, ps0
    use shr_vmath_mod,  only: shr_vmath_log
    use phys_gmean,     only: gmean
    !---------------------------------------------------------------------------
    implicit none
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer  :: lchnk
    real(r8) :: qbot                 ! bottom level q before change
    real(r8) :: qbotm1               ! bottom-1 level q before change
    real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
    real(r8) :: qmavl                ! available q at level pver-1
    real(r8) :: ke(pcols,begchunk:endchunk)   
    real(r8) :: se(pcols,begchunk:endchunk)   
    real(r8) :: ke_glob(1),se_glob(1)
    real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer
    integer  :: m, i, k, ncol
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)
    !---------------------------------------------------------------------------

    ! Evaluate derived quantities
    !$omp parallel do private (lchnk, ncol, k, i, zvirv, pbuf_chnk)
    do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      do k = 1,nlev
        do i = 1,ncol
          phys_state(lchnk)%pint(i,k)=hyai(k)*ps0+hybi(k)*phys_state(lchnk)%ps(i)
          phys_state(lchnk)%pmid(i,k)=hyam(k)*ps0+hybm(k)*phys_state(lchnk)%ps(i)
        end do
        call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k), &
                           phys_state(lchnk)%lnpint(1:ncol,k),ncol)
        call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k), &
                           phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
      end do
      do i = 1,ncol
        phys_state(lchnk)%pint(i,pverp)=hyai(pverp)*ps0+hybi(pverp)*phys_state(lchnk)%ps(i)
      end do
      call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp), &
                         phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)

      do k = 1,nlev
        do i = 1,ncol
          phys_state(lchnk)%pdel (i,k)  = phys_state(lchnk)%pint(i,k+1) &
                                         -phys_state(lchnk)%pint(i,k)
          phys_state(lchnk)%rpdel(i,k)  = 1._r8/phys_state(lchnk)%pdel(i,k)
          phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                                          /phys_state(lchnk)%pmid(i,k))**cappa
        end do
      end do

      !----------------------------------------------------
      ! Need to fill zvirv 2D variables to be 
      ! compatible with geopotential_t interface
      !----------------------------------------------------
      zvirv(:,:) = zvir

      ! Compute initial geopotential heights
      call geopotential_t(phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  ,&
                          phys_state(lchnk)%pint  , phys_state(lchnk)%pmid    ,&
                          phys_state(lchnk)%pdel  , phys_state(lchnk)%rpdel   ,&
                          phys_state(lchnk)%t     , phys_state(lchnk)%q(:,:,1),&
                          rairv(:,:,lchnk)        , gravit, zvirv             ,&
                          phys_state(lchnk)%zi    , phys_state(lchnk)%zm      ,&
                          ncol)
          
       ! Compute initial dry static energy s = g*z + c_p*T
       do k = 1, pver
          do i=1,ncol
             phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k)    &
                                      + gravit*phys_state(lchnk)%zm(i,k)  &
                                      + phys_state(lchnk)%phis(i)
          end do
       end do

       ! NOTE:  if a tracer is marked "dry", that means physics wants it dry
       !        if dycore advects it wet, it should be converted here 
       !        FV dycore does this, and in physics/cam/tphysac.F90 it will
       !        be converted back to wet, BUT ONLY FOR FV dycore
       !
       !        EUL: advects all tracers (except q1) as dry.  so it never 
       !        calls this.
       !
       !        SE:  we should follow FV and advect all tracers wet (especially
       !        since we will be switching to conservation form of advection).  
       !        So this is broken since dry tracers will never get converted 
       !        back to wet. (in APE, all tracers are wet, so it is ok for now)  
       !
       ! Convert dry type constituents from moist to dry mixing ratio
       call set_state_pdry(phys_state(lchnk))	! First get dry pressure to use for this timestep
       call set_wet_to_dry(phys_state(lchnk)) ! Dynamics had moist, physics wants dry.

       ! Compute energy and water integrals of input state
       pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
       call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)

	
#if 0
       ke(:,lchnk) = 0._r8
       se(:,lchnk) = 0._r8
      ! wv = 0._r8
      ! wl = 0._r8
      ! wi = 0._r8
       do k = 1, pver
          do i = 1, ncol
             ke(i,lchnk) = ke(i,lchnk) + ( 0.5_r8*(phys_state(lchnk)%u(i,k)**2 + &
                  phys_state(lchnk)%v(i,k)**2)*phys_state(lchnk)%pdel(i,k) )/gravit
             se(i,lchnk) = se(i,lchnk) + phys_state(lchnk)%s(i,k         )*phys_state(lchnk)%pdel(i,k)/gravit
            ! wv = wv + phys_state(lchnk)%q(i,k,1       )*phys_state(lchnk)%pdel(i,k)
            ! wl = wl + phys_state(lchnk)%q(i,k,ixcldliq)*phys_state(lchnk)%pdel(i,k)
            ! wi = wi + phys_state(lchnk)%q(i,k,ixcldice)*phys_state(lchnk)%pdel(i,k)
          end do
       end do
#endif 
    end do ! lchnk

#if 0
    ! This wont match SE exactly.  SE computes KE at half levels
    ! SE includes cp_star (physics SE uses cp )
    ! CAM stdout of total energy also includes latent energy of Q,Q1,Q2
    ! making it a little larger
    call gmean(ke,ke_glob,1)
    call gmean(se,se_glob,1)
    if (masterproc) then
       write(iulog,'(a,e20.8,a,e20.8)') 'KE = ',ke_glob(1),' SE = ',se_glob(1)
       write(iulog,'(a,e20.8)') 'TOTE = ',ke_glob(1)+se_glob(1)
    end if ! masterproc
#endif

  end subroutine derived_phys
  !=================================================================================================
  !=================================================================================================
end module dp_coupling
