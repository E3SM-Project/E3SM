!
!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------
module dp_coupling
  use constituents,   only: pcnst, cnst_name
  use cam_history,    only: outfld, write_inithist, hist_fld_active
  use dimensions_mod, only: np, npsq, nelemd, nlev
  use dof_mod,        only: UniquePoints, PutUniquePoints
  use dyn_comp,       only: dyn_export_t, dyn_import_t, TimeLevel
  use dyn_grid,       only: get_gcol_block_d
  use element_mod,    only: element_t
  use kinds,          only: real_kind, int_kind
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use physics_types,  only: physics_state, physics_tend
  
  use phys_grid,      only: get_ncols_p, get_gcol_all_p, block_to_chunk_send_pters, transpose_block_to_chunk, &
       block_to_chunk_recv_pters, chunk_to_block_send_pters, transpose_chunk_to_block, chunk_to_block_recv_pters
  use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
  use element_mod,    only: element_t
  use control_mod,    only: smooth_phis_numcycle
  use cam_logfile,    only : iulog
  use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
  use spmd_utils,   only: mpicom, iam
  use perf_mod,    only : t_startf, t_stopf, t_barrierf
  use parallel_mod, only : par
  private
  public :: d_p_coupling, p_d_coupling
!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
  subroutine d_p_coupling(phys_state, phys_tend,  pbuf2d, dyn_out)
    use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk, &
                              pbuf_get_field
    use shr_vmath_mod,  only: shr_vmath_exp
    use time_manager,   only: is_first_step
    use viscosity_mod,  only: compute_zeta_C0
    use cam_abortutils,     only: endrun
    use gravity_waves_sources, only: gws_src_fnct
    use dyn_comp,       only: frontgf_idx, frontga_idx
    use phys_control,   only: use_gw_front
    implicit none
!-----------------------------------------------------------------------
! !INPUT PARAMETERS:
!
    type(dyn_export_t), intent(inout) :: dyn_out    ! dynamics export 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
! !OUTPUT PARAMETERS:

    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    

! LOCAL VARIABLES
    type(element_t), pointer :: elem(:)               ! pointer to dyn_out element array
    integer (kind=int_kind)  :: ie               ! indices over elements
    integer (kind=int_kind)  :: lchnk, icol, ilyr      ! indices over chunks, columns, layers
    real (kind=real_kind)    :: ps_tmp(npsq,nelemd)         ! temporary array to hold ps
    real (kind=real_kind)    :: phis_tmp(npsq,nelemd)       ! temporary array to hold phis  
    real (kind=real_kind)    :: T_tmp(npsq,pver,nelemd)     ! temporary array to hold T
    real (kind=real_kind)    :: uv_tmp(npsq,2,pver,nelemd)     ! temporary array to hold u and v
    real (kind=real_kind)    :: q_tmp(npsq,pver,pcnst,nelemd) ! temporary to hold advected constituents
    real (kind=real_kind)    :: omega_tmp(npsq,pver,nelemd) ! temporary array to hold omega

    ! Frontogenesis
    real (kind=real_kind), allocatable :: frontgf(:,:,:) ! temporary arrays to hold frontogenesis
    real (kind=real_kind), allocatable :: frontga(:,:,:) !   function (frontgf) and angle (frontga)
    ! Pointers to pbuf
    real (kind=r8),        pointer     :: pbuf_frontgf(:,:)
    real (kind=r8),        pointer     :: pbuf_frontga(:,:)

    integer :: ncols,i,j,ierr

    integer (kind=int_kind)  :: ioff, m
    integer :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: bpter(npsq,0:pver)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

    real (kind=real_kind), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers
    integer :: tl_f

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    !----------------------------------------------------------------------

    nullify(pbuf_chnk)
    nullify(pbuf_frontgf)
    nullify(pbuf_frontga)

    if (use_gw_front) then

       allocate(frontgf(npsq,pver,nelemd), stat=ierr)
       if (ierr /= 0) call endrun("dp_coupling: Allocate of frontgf failed.")

       allocate(frontga(npsq,pver,nelemd), stat=ierr)
       if (ierr /= 0) call endrun("dp_coupling: Allocate of frontga failed.")

    end if

    if( iam < par%nprocs) then

       elem => dyn_out%elem

       tl_f = TimeLevel%n0  ! time split physics (with forward-in-time RK)

       call t_startf('UniquePoints')
       do ie=1,nelemd
          ncols = elem(ie)%idxP%NumUniquePts
          call UniquePoints(elem(ie)%idxP, elem(ie)%state%ps_v(:,:,tl_f), ps_tmp(1:ncols,ie))
          call UniquePoints(elem(ie)%idxP, nlev, elem(ie)%state%T(:,:,:,tl_f), T_tmp(1:ncols,:,ie))
          call UniquePoints(elem(ie)%idxP, 2, nlev, elem(ie)%state%V(:,:,:,:,tl_f), uv_tmp(1:ncols,:,:,ie))
          call UniquePoints(elem(ie)%idxP, nlev, elem(ie)%derived%omega_p, omega_tmp(1:ncols,:,ie))

          call UniquePoints(elem(ie)%idxP, elem(ie)%state%phis, phis_tmp(1:ncols,ie))
          call UniquePoints(elem(ie)%idxP, nlev,pcnst, elem(ie)%state%Q(:,:,:,:), Q_tmp(1:ncols,:,:,ie))
       end do
       call t_stopf('UniquePoints')

       if (use_gw_front) call gws_src_fnct(elem, tl_f, frontgf, frontga)
    else
       ps_tmp(:,:) = 0._r8
       T_tmp(:,:,:) = 0._r8
       uv_tmp(:,:,:,:) = 0._r8
       omega_tmp(:,:,:) = 0._r8
       if (use_gw_front) then
          frontgf(:,:,:) = 0._r8
          frontga(:,:,:) = 0._r8
       end if
       phis_tmp(:,:) = 0._r8
       Q_tmp(:,:,:,:) = 0._r8
    endif !! iam .lt. par%nprocs

    call t_startf('dpcopy')
    if (local_dp_map) then

!$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, ilyr, m, pbuf_chnk, pbuf_frontgf, pbuf_frontga)
       do lchnk=begchunk,endchunk
          ncols=get_ncols_p(lchnk)
          call get_gcol_all_p(lchnk,pcols,pgcols)

          pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

          if (use_gw_front) then
             call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
             call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
          end if

          do icol=1,ncols
             call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
             ie = idmb3(1)
             ioff=idmb2(1)

             phys_state(lchnk)%ps(icol)=ps_tmp(ioff,ie)
             phys_state(lchnk)%phis(icol)=phis_tmp(ioff,ie)
             do ilyr=1,pver
                phys_state(lchnk)%t(icol,ilyr)=T_tmp(ioff,ilyr,ie)	   
                phys_state(lchnk)%u(icol,ilyr)=uv_tmp(ioff,1,ilyr,ie)
                phys_state(lchnk)%v(icol,ilyr)=uv_tmp(ioff,2,ilyr,ie)
                phys_state(lchnk)%omega(icol,ilyr)=omega_tmp(ioff,ilyr,ie)

                if (use_gw_front) then
                   pbuf_frontgf(icol,ilyr) = frontgf(ioff,ilyr,ie)
                   pbuf_frontga(icol,ilyr) = frontga(ioff,ilyr,ie)
                endif
             end do

             do m=1,pcnst
                do ilyr=1,pver
                   phys_state(lchnk)%q(icol,ilyr,m)=Q_tmp(ioff,ilyr,m,ie)
                end do
             end do
          end do
       end do

    else  ! .not. local_dp_map

       tsize = 4 + pcnst
       if (use_gw_front) tsize = tsize + 2

       allocate(bbuffer(tsize*block_buf_nrecs))
       allocate(cbuffer(tsize*chunk_buf_nrecs))
       if(iam .lt. par%nprocs) then
!$omp parallel do private (ie, bpter, icol, ilyr, m)
          do ie=1,nelemd

             call block_to_chunk_send_pters(elem(ie)%GlobalID,npsq,pver+1,tsize,bpter)

             do icol=1,elem(ie)%idxP%NumUniquePts
                
                bbuffer(bpter(icol,0)+2:bpter(icol,0)+tsize-1) = 0.0_r8
                
                bbuffer(bpter(icol,0))   = ps_tmp(icol,ie)
                bbuffer(bpter(icol,0)+1) = phis_tmp(icol,ie)

                do ilyr=1,pver
                   bbuffer(bpter(icol,ilyr))   = T_tmp(icol,ilyr,ie)
                   bbuffer(bpter(icol,ilyr)+1) = uv_tmp(icol,1,ilyr,ie)
                   bbuffer(bpter(icol,ilyr)+2) = uv_tmp(icol,2,ilyr,ie)
                   bbuffer(bpter(icol,ilyr)+3) = omega_tmp(icol,ilyr,ie)

                   if (use_gw_front) then
                      bbuffer(bpter(icol,ilyr)+4) = frontgf(icol,ilyr,ie)
                      bbuffer(bpter(icol,ilyr)+5) = frontga(icol,ilyr,ie)
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

!$omp parallel do private (lchnk, ncols, cpter, icol, ilyr, m, pbuf_chnk, pbuf_frontgf, pbuf_frontga)
       do lchnk = begchunk,endchunk
          ncols = phys_state(lchnk)%ncol

          pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

          if (use_gw_front) then
             call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
             call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
          end if

          call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

          do icol=1,ncols

             phys_state(lchnk)%ps   (icol)     = cbuffer(cpter(icol,0))
             phys_state(lchnk)%phis (icol)     = cbuffer(cpter(icol,0)+1)

             do ilyr=1,pver

                phys_state(lchnk)%t     (icol,ilyr)   = cbuffer(cpter(icol,ilyr))
                phys_state(lchnk)%u     (icol,ilyr)   = cbuffer(cpter(icol,ilyr)+1)
                phys_state(lchnk)%v     (icol,ilyr)   = cbuffer(cpter(icol,ilyr)+2)
                phys_state(lchnk)%omega (icol,ilyr)   = cbuffer(cpter(icol,ilyr)+3)

                if (use_gw_front) then
                   pbuf_frontgf(icol,ilyr) = cbuffer(cpter(icol,ilyr)+4)
                   pbuf_frontga(icol,ilyr) = cbuffer(cpter(icol,ilyr)+5)
                endif

                do m=1,pcnst
                   phys_state(lchnk)%q  (icol,ilyr,m) = cbuffer(cpter(icol,ilyr)+tsize-pcnst-1+m)
                end do

             end do

          end do

       end do

       deallocate( bbuffer )
       deallocate( cbuffer )

    end if
    call t_stopf('dpcopy')

    call t_startf('derived_phys')
    call derived_phys(phys_state,phys_tend,pbuf2d)
    call t_stopf('derived_phys')

!$omp parallel do private (lchnk, ncols, ilyr, icol)
    do lchnk=begchunk,endchunk
       ncols=get_ncols_p(lchnk)
       do ilyr=1,pver
          do icol=1,ncols
             phys_state(lchnk)%omega(icol,ilyr)=phys_state(lchnk)%omega(icol,ilyr)*phys_state(lchnk)%pmid(icol,ilyr)
          end do
       end do
    end do


   if (write_inithist() ) then
      do lchnk=begchunk,endchunk
         call outfld('T&IC',phys_state(lchnk)%t,pcols,lchnk)
         call outfld('U&IC',phys_state(lchnk)%u,pcols,lchnk)
         call outfld('V&IC',phys_state(lchnk)%v,pcols,lchnk)
         call outfld('PS&IC',phys_state(lchnk)%ps,pcols,lchnk)
         do m=1,pcnst
            call outfld(trim(cnst_name(m))//'&IC',phys_state(lchnk)%q(1,1,m), pcols,lchnk)
         end do
      end do
   endif
   
       
  end subroutine d_p_coupling

  subroutine p_d_coupling(phys_state, phys_tend,  dyn_in)
    use shr_vmath_mod, only: shr_vmath_log
    use cam_control_mod, only : adiabatic
    implicit none

! !INPUT PARAMETERS:
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend),  intent(inout), dimension(begchunk:endchunk) :: phys_tend
    

! !OUTPUT PARAMETERS:
    type(dyn_import_t),  intent(inout)   :: dyn_in

    ! LOCAL VARIABLES
    integer :: ic , ncols                            ! index
    type(element_t), pointer :: elem(:)               ! pointer to dyn_in element array
    integer (kind=int_kind)  :: ie, iep               ! indices over elements
    integer (kind=int_kind)  :: lchnk, icol, ilyr      ! indices over chunks, columns, layers

    real (kind=real_kind)    :: T_tmp(npsq,pver,nelemd)       ! temporary array to hold T
    real (kind=real_kind)    :: uv_tmp(npsq,2,pver,nelemd)    ! temporary array to hold uv
!   real (kind=real_kind)    :: omega_tmp(npsq,pver,nelemd)   ! temporary array to hold omega
    real (kind=real_kind)    :: q_tmp(npsq,pver,pcnst,nelemd) ! temporary array to hold q
    integer (kind=int_kind)  :: ioff, m, i, j, k
    integer(kind=int_kind)   :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)

    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for packing data
    integer :: bpter(npsq,0:pver)    ! offsets into block buffer for unpacking data

    real (kind=real_kind), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

    if (iam .lt. par%nprocs) then
       elem => dyn_in%elem
    else
       nullify(elem)
    end if

    T_tmp=0.0_r8
    uv_tmp=0.0_r8
    q_tmp=0.0_r8

    if(adiabatic) return

    call t_startf('pd_copy')
    if(local_dp_map) then

!$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, ilyr, m)
       do lchnk=begchunk,endchunk
          ncols=get_ncols_p(lchnk)
          call get_gcol_all_p(lchnk,pcols,pgcols)

          do icol=1,ncols
             call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
             ie = idmb3(1)
             ioff=idmb2(1)

             do ilyr=1,pver
                T_tmp(ioff,ilyr,ie)      = phys_tend(lchnk)%dtdt(icol,ilyr)
                uv_tmp(ioff,1,ilyr,ie)   = phys_tend(lchnk)%dudt(icol,ilyr)
                uv_tmp(ioff,2,ilyr,ie)   = phys_tend(lchnk)%dvdt(icol,ilyr)

                do m=1,pcnst
                   q_tmp(ioff,ilyr,m,ie) = phys_state(lchnk)%q(icol,ilyr,m)
                end do
             end do

       	  end do	

       end do

    else

       tsize = 3 + pcnst

       allocate( bbuffer(tsize*block_buf_nrecs) )
       allocate( cbuffer(tsize*chunk_buf_nrecs) )

!$omp parallel do private (lchnk, ncols, cpter, i, icol, ilyr, m)
       do lchnk = begchunk,endchunk
          ncols = get_ncols_p(lchnk)

          call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)

          do i=1,ncols
             cbuffer(cpter(i,0):cpter(i,0)+2+pcnst) = 0.0_r8
          end do

          do icol=1,ncols

             do ilyr=1,pver
                cbuffer   (cpter(icol,ilyr))     = phys_tend(lchnk)%dtdt(icol,ilyr)
                cbuffer   (cpter(icol,ilyr)+1)   = phys_tend(lchnk)%dudt(icol,ilyr)
                cbuffer   (cpter(icol,ilyr)+2)   = phys_tend(lchnk)%dvdt(icol,ilyr)

                do m=1,pcnst
                   cbuffer(cpter(icol,ilyr)+2+m) = phys_state(lchnk)%q(icol,ilyr,m)
                end do
             end do

          end do

       end do

       call t_barrierf('sync_chk_to_blk', mpicom)
       call t_startf ('chunk_to_block')
       call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
       call t_stopf  ('chunk_to_block')
       if(iam < par%nprocs) then
!$omp parallel do private (ie, bpter, icol, ilyr, m)
          do ie=1,nelemd

             call chunk_to_block_recv_pters(elem(ie)%GlobalID,npsq,pver+1,tsize,bpter)

             do icol=1,elem(ie)%idxP%NumUniquePts

                do ilyr=1,pver

                   T_tmp   (icol,ilyr,ie)   = bbuffer(bpter(icol,ilyr))
                   uv_tmp  (icol,1,ilyr,ie) = bbuffer(bpter(icol,ilyr)+1)
                   uv_tmp  (icol,2,ilyr,ie) = bbuffer(bpter(icol,ilyr)+2)

                   do m=1,pcnst
                      q_tmp(icol,ilyr,m,ie) = bbuffer(bpter(icol,ilyr)+2+m)
                   end do

                end do
                
             end do
             
          end do
       endif
       deallocate( bbuffer )
       deallocate( cbuffer )
       
    end if
    call t_stopf('pd_copy')
    if(iam < par%nprocs) then
       call t_startf('putUniquePoints')
       do ie=1,nelemd
          ncols = elem(ie)%idxP%NumUniquePts
          call putUniquePoints(elem(ie)%idxP, nlev, T_tmp(1:ncols,:,ie), elem(ie)%derived%fT(:,:,:,1))
          call putUniquePoints(elem(ie)%idxP, 2, nlev, uv_tmp(1:ncols,:,:,ie), &
               elem(ie)%derived%fM(:,:,:,:,1))
          call putUniquePoints(elem(ie)%idxP, nlev,pcnst, q_tmp(1:ncols,:,:,ie), &
               elem(ie)%derived%fQ(:,:,:,:,1))
       end do
       call t_stopf('putUniquePoints')
    end if
  end subroutine p_d_coupling

  subroutine derived_phys(phys_state, phys_tend, pbuf2d)
    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk
    use constituents,  only: qmin
    use physconst,     only: cpair, gravit, rair, zvir, cappa, rairv
    use spmd_utils,    only: masterproc
    use ppgrid,        only: pver
    use geopotential,  only: geopotential_t
    use physics_types, only: set_state_pdry, set_wet_to_dry
    use check_energy,  only: check_energy_timestep_init
    use hycoef,   only : hyam, hybm, hyai, hybi, ps0
    use shr_vmath_mod, only: shr_vmath_log
    use phys_gmean,      only: gmean


    implicit none
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(physics_buffer_desc),      pointer     :: pbuf2d(:,:)

    integer :: lchnk
    real(r8) :: qbot                 ! bottom level q before change
    real(r8) :: qbotm1               ! bottom-1 level q before change
    real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
    real(r8) :: qmavl                ! available q at level pver-1

    real(r8) :: ke(pcols,begchunk:endchunk)   
    real(r8) :: se(pcols,begchunk:endchunk)   
    real(r8) :: ke_glob(1),se_glob(1)
    real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer

    integer :: m, i, k, ncol

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

!
! Evaluate derived quantities
!
!$omp parallel do private (lchnk, ncol, k, i, zvirv, pbuf_chnk)
    do lchnk = begchunk,endchunk
       ncol = get_ncols_p(lchnk)
       do k=1,nlev
          do i=1,ncol
             phys_state(lchnk)%pint(i,k)=hyai(k)*ps0+hybi(k)*phys_state(lchnk)%ps(i)
             phys_state(lchnk)%pmid(i,k)=hyam(k)*ps0+hybm(k)*phys_state(lchnk)%ps(i)
          end do
          call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k),phys_state(lchnk)%lnpint(1:ncol,k),ncol)
          call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k),phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
       end do
       do i=1,ncol
          phys_state(lchnk)%pint(i,pverp)=hyai(pverp)*ps0+hybi(pverp)*phys_state(lchnk)%ps(i)
       end do
       call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp),phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)


       do k=1,nlev
          do i=1,ncol
             phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pint(i,k+1) - phys_state(lchnk)%pint(i,k)
             phys_state(lchnk)%rpdel(i,k) = 1._r8/phys_state(lchnk)%pdel(i,k)
             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                                             / phys_state(lchnk)%pmid(i,k))**cappa
          end do
       end do

!-----------------------------------------------------------------------------------
!  Need to fill zvirv 2D variables to be compatible with geopotential_t interface
!-----------------------------------------------------------------------------------
       zvirv(:,:) = zvir
!
! Compute initial geopotential heights
!

       call geopotential_t (phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
            phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
            phys_state(lchnk)%t     , phys_state(lchnk)%q(:,:,1), rairv(:,:,lchnk),  gravit,  zvirv       , &
            phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol                )
          
! Compute initial dry static energy, include surface geopotential
       do k = 1, pver
          do i=1,ncol
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
       call set_state_pdry(phys_state(lchnk))	 ! First get dry pressure to use for this timestep
       call set_wet_to_dry(phys_state(lchnk))    ! Dynamics had moist, physics wants dry.

! Compute energy and water integrals of input state
       pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
       call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)

	
#if 0
       ke(:,lchnk) = 0._r8
       se(:,lchnk) = 0._r8
!       wv = 0._r8
!       wl = 0._r8
!       wi = 0._r8
       do k = 1, pver
          do i = 1, ncol
             ke(i,lchnk) = ke(i,lchnk) + ( 0.5_r8*(phys_state(lchnk)%u(i,k)**2 + &
                  phys_state(lchnk)%v(i,k)**2)*phys_state(lchnk)%pdel(i,k) )/gravit
             se(i,lchnk) = se(i,lchnk) + phys_state(lchnk)%s(i,k         )*phys_state(lchnk)%pdel(i,k)/gravit
!             wv = wv + phys_state(lchnk)%q(i,k,1       )*phys_state(lchnk)%pdel(i,k)
!             wl = wl + phys_state(lchnk)%q(i,k,ixcldliq)*phys_state(lchnk)%pdel(i,k)
!             wi = wi + phys_state(lchnk)%q(i,k,ixcldice)*phys_state(lchnk)%pdel(i,k)
          end do
       end do
#endif 
    end do

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
    endif
#endif

  end subroutine derived_phys


end module dp_coupling
