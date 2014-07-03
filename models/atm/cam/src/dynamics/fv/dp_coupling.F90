module dp_coupling
!BOP
!
! !MODULE: dp_coupling --- dynamics-physics coupling module
!
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use rgrid,         only: nlon
   use ppgrid,        only: pcols, pver, pverp
   use phys_grid
   
   use physics_types, only: physics_state, physics_tend
   use constituents,  only: pcnst, qmin
   use physconst,     only: cpair, gravit, rair, zvir, cpairv, rairv
   use geopotential,  only: geopotential_t
   use check_energy,  only: check_energy_timestep_init
   use dynamics_vars, only: T_FVDYCORE_GRID
   use dyn_comp,      only: dyn_import_t, dyn_export_t
   use abortutils,    only: endrun
#if defined ( SPMD )
   use spmd_dyn,      only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
#endif
   use perf_mod
   use cam_logfile,   only: iulog
   use phys_control,  only: do_waccm_phys, waccmx_is !WACCM query functions

!--------------------------------------------
!  Variables needed for WACCM-X
!--------------------------------------------
   use constituents,  only: cnst_get_ind, cnst_mw  !Needed to access constituent molecular weights
   use shr_const_mod, only: shr_const_rgas         !Gas constant
!
! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC d_p_coupling, p_d_coupling

!
! !DESCRIPTION:
!
!      This module provides 
!
!      \begin{tabular}{|l|l|} \hline \hline
!        d\_p\_coupling    &  dynamics output to physics input \\ \hline
!        p\_d\_coupling    &  physics output to dynamics input \\ \hline 
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   00.06.01   Boville    Creation
!   01.10.01   Lin        Various revisions
!   01.03.26   Sawyer     Added ProTeX documentation
!   01.06.27   Mirin      Separate noncoupling coding into new routines
!   01.07.13   Mirin      Some support for multi-2D decompositions
!   02.03.01   Worley     Support for nontrivial physics remapping
!   03.03.28   Boville    set all physics_state elements, add check_energy_timestep_init
!   03.08.13   Sawyer     Removed ghost N1 region in u3sxy
!   05.06.28   Sawyer     Simplified interfaces -- only XY decomposition 
!   05.10.25   Sawyer     Extensive refactoring, dyn_interface
!   05.11.10   Sawyer     Now using dyn_import/export_t containers
!   06.07.01   Sawyer     Transitioned constituents to T_TRACERS
!
!EOP
!-----------------------------------------------------------------------

   private 
   real(r8), parameter ::  D0_5                    =  0.5_r8
   real(r8), parameter ::  D1_0                    =  1.0_r8

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: d_p_coupling --- convert dynamics output to physics input
!
! !INTERFACE: 
  subroutine d_p_coupling(grid, phys_state, phys_tend,  pbuf2d, dyn_out)

! !USES:
    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk
    use constituents,  only: cnst_get_type_byind
    use physics_types, only: set_state_pdry, set_wet_to_dry, &
         physics_state_dycore_alloc

    use pmgrid, only : plev, plevp
    use ctem, only   : ctem_diags, do_circulation_diags
    use gravity_waves_sources, only: gws_src_fnct
    use physconst,    only: physconst_update
    use shr_const_mod, only: shr_const_rwv

!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
! !INPUT PARAMETERS:
!
    type(T_FVDYCORE_GRID), intent(in) :: grid ! FV Dynamics grid
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(dyn_export_t), intent(in)    :: dyn_out  ! dynamics export 

! !OUTPUT PARAMETERS:

    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    

! !DESCRIPTION:
!
!   Coupler for converting dynamics output variables into physics 
!   input variables
!
! !REVISION HISTORY:
!   00.06.01   Boville    Creation
!   01.07.13   AAM        Some support for multi-2D decompositions
!   02.03.01   Worley     Support for nontrivial physics remapping
!   02.05.02   Sawyer     u3s made inout due to ghosting in d2a3dikj
!   03.08.05   Sawyer     Removed pe11k, pe11kln (for defunct Rayl fric)
!   04.08.29   Eaton      Added lat, lon coords to physics_state type
!   05.06.28   Sawyer     Simplified interface -- on XY decomp vars.
!   05.07.06   Sawyer     Added dyn_state as argument
!   05.10.31   Sawyer     Refactoring, replaced dyn_state by dyn_interface
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

! Variables from dynamics export container
    real(r8), pointer :: phisxy(:,:)              ! surface geopotential
    real(r8), pointer :: psxy (:,:)               ! surface pressure
    real(r8), pointer :: u3sxy(:,:,:)             ! u-wind on d-grid
    real(r8), pointer :: v3sxy(:,:,:)             ! v-wind on d-grid
    real(r8), pointer :: ptxy (:,:,:)             ! Virtual pot temp
    real(r8), pointer :: tracer(:,:,:,:)          ! constituents
    real(r8), pointer :: omgaxy(:,:,:)            ! vertical velocity
    real(r8), pointer :: pexy  (:,:,:)            ! edge pressure
    real(r8), pointer :: pelnxy(:,:,:)            ! log(pe)
    real(r8), pointer :: pkxy  (:,:,:)            ! pe**cappa
    real(r8), pointer :: pkzxy (:,:,:)            ! f-v mean of pk

    integer :: i,ib,j,k,m,lchnk      ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: blksiz                ! number of columns in 2D block
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer, allocatable, dimension(:,:) :: bpter
                                     ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

    real(r8) :: rlat(pcols)          ! array of latitudes (radians)
    real(r8) :: rlon(pcols)          ! array of longitudes (radians)
    real(r8) :: qmavl                ! available q at level pver-1
    real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
    real(r8) :: qbot                 ! bottom level q before change
    real(r8) :: qbotm1               ! bottom-1 level q before change
    real(r8) :: pic(pcols)           ! ps**cappa
    real(r8) :: fraction
    real(r8), allocatable :: u3(:, :, :)       ! u-wind on a-grid
    real(r8), allocatable :: v3(:, :, :)       ! v-wind on a-grid
    real(r8), allocatable, dimension(:) :: bbuffer, cbuffer
                                     ! transpose buffers

    real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer

    integer  :: im, jm, km, kmp1, iam
    integer  :: ifirstxy, ilastxy, jfirstxy, jlastxy
    integer  :: ic, jc
    integer  :: astat
    integer  :: boff
    logical, save :: debug_adjust_print = .true. ! true => print out tracer adjustment msgs

!--------------------------------------------
!  Variables needed for WACCM
!--------------------------------------------

    ! frontogenesis function for gravity wave drag
    real(r8), allocatable :: frontgf(:,:,:)
    ! frontogenesis angle for gravity wave drag
    real(r8), allocatable :: frontga(:,:,:)
    ! needed for qbo
    real(r8) :: uzm(plev,grid%jfirstxy:grid%jlastxy)

!--------------------------------------------
!  Variables needed for WACCM-X
!--------------------------------------------
    integer  :: ixo, ixo2, ixh, ixh2, ixn  ! indices into state structure for O, O2, H, H2, and N
    real(r8) :: mmrSum_O_O2_H                ! Sum of mass mixing ratios for O, O2, and H
    real(r8), parameter :: mmrMin=1.e-20_r8  ! lower limit of o2, o, and h mixing ratios
    real(r8), parameter :: N2mmrMin=1.e-6_r8 ! lower limit of o2, o, and h mixing ratios

#if (! defined SPMD)
    integer  :: block_buf_nrecs = 0
    integer  :: chunk_buf_nrecs = 0
    logical  :: local_dp_map=.true. 
#endif
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

!---------------------------End Local workspace-------------------------

    if (do_waccm_phys()) then
       ! Currently, only WACCM uses dycore-initialized fields in physics_state.

!$omp parallel do private (lchnk)
       do lchnk = begchunk, endchunk
          if (.not. phys_state(lchnk)%dycore_alloc) &
               call physics_state_dycore_alloc(phys_state(lchnk))
       end do

       allocate(frontgf(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy), stat=astat)
       if( astat /= 0 ) then
          write(iulog,*) 'd_p_coupling: failed to allocate frontgf; error = ',astat
          call endrun
       end if
       allocate(frontga(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy), stat=astat)
       if( astat /= 0 ) then
          write(iulog,*) 'd_p_coupling: failed to allocate frontga; error = ',astat
          call endrun
       end if
    end if

    fraction = 0.1_r8

    phisxy   => dyn_out%phis
    psxy     => dyn_out%ps
    u3sxy    => dyn_out%u3s
    v3sxy    => dyn_out%v3s
    ptxy     => dyn_out%pt
    tracer   => dyn_out%tracer

    omgaxy   => dyn_out%omga
    pexy     => dyn_out%pe
    pelnxy   => dyn_out%peln
    pkxy     => dyn_out%pk
    pkzxy    => dyn_out%pkz

    im       = grid%im
    jm       = grid%jm
    km       = grid%km
    kmp1     = km + 1

    ifirstxy = grid%ifirstxy
    ilastxy  = grid%ilastxy
    jfirstxy = grid%jfirstxy
    jlastxy  = grid%jlastxy

    iam      = grid%iam
!-----------------------------------------------------------------------
! Transform dynamics staggered winds to physics grid (D=>A)
!-----------------------------------------------------------------------

    call t_startf ('d2a3dikj')
    allocate (u3(ifirstxy:ilastxy, km, jfirstxy:jlastxy))
    allocate (v3(ifirstxy:ilastxy, km, jfirstxy:jlastxy))

    if (iam .lt. grid%npes_xy) then
       call d2a3dikj( grid, u3sxy,  v3sxy, u3, v3 )
    end if  ! (iam .lt. grid%npes_xy)

    call t_stopf  ('d2a3dikj')

    if ( do_circulation_diags ) then
       call t_startf('DP_CPLN_ctem')
       call ctem_diags( u3, v3, omgaxy, ptxy(:,jfirstxy:jlastxy,:), tracer(:,jfirstxy:jlastxy,:,1), &
                         psxy, pexy, grid, uzm )
       call t_stopf('DP_CPLN_ctem')
    endif
    if (do_waccm_phys()) then
       call t_startf('DP_CPLN_gw_sources')
       call gws_src_fnct (u3,v3,ptxy,  tracer(:,jfirstxy:jlastxy,:,1), pexy, grid, frontgf, frontga)
       call t_stopf('DP_CPLN_gw_sources')
    end if
         
!-----------------------------------------------------------------------
! Copy data from dynamics data structure to physics data structure
!-----------------------------------------------------------------------
has_local_map : &
    if (local_dp_map) then

!$omp parallel do private (lchnk, ncol, i, k, m, ic, jc, lons, lats, pic)
chnk_loop1 : &
       do lchnk = begchunk,endchunk
          ncol = phys_state(lchnk)%ncol
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)

          do i=1,ncol
             ic = lons(i)
             jc = lats(i)
             phys_state(lchnk)%ps(i)   = psxy(ic,jc)
             phys_state(lchnk)%phis(i) = phisxy(ic,jc)
             pic(i) = pkxy(ic,jc,pver+1)
          enddo
          do k=1,km
             do i=1,ncol
                ic = lons(i)
                jc = lats(i)
                phys_state(lchnk)%u    (i,k) = u3(ic,k,jc)
                phys_state(lchnk)%v    (i,k) = v3(ic,k,jc)
                phys_state(lchnk)%omega(i,k) = omgaxy(ic,k,jc)
                if (do_waccm_phys()) then
                   phys_state(lchnk)%frontgf(i,k) = frontgf(ic,k,jc)
                   phys_state(lchnk)%frontga(i,k) = frontga(ic,k,jc)
                   phys_state(lchnk)%uzm(i,k)     = uzm(k,jc)
                endif
                phys_state(lchnk)%t    (i,k) = ptxy(ic,jc,k) / (D1_0 + zvir*tracer(ic,jc,k,1))
                phys_state(lchnk)%exner(i,k) = pic(i) / pkzxy(ic,jc,k) 
             end do
          end do

          do k=1,kmp1
             do i=1,ncol
!
! edge-level pressure arrays: copy from the arrays computed by dynpkg
!
                ic = lons(i)
                jc = lats(i)
                phys_state(lchnk)%pint  (i,k) = pexy  (ic,k,jc)
                phys_state(lchnk)%lnpint(i,k) = pelnxy(ic,k,jc)
             end do
          end do

!
! Copy constituents
! Dry types converted from moist to dry m.r. at bottom of this routine
!
          do m=1,pcnst
             do k=1,km
                do i=1,ncol
                   phys_state(lchnk)%q(i,k,m) = &
                      tracer(lons(i),lats(i),k,m)
                end do
             end do
          end do
 
       end do chnk_loop1

    else has_local_map

       tsize = 7 + pcnst
       boff  = 6
       if (do_waccm_phys()) then
          tsize = tsize+3
          boff  = boff+3
       end if
 
       blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
       allocate( bpter(blksiz,0:km),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'd_p_coupling: failed to allocate bpter; error = ',astat
          call endrun
       end if
       allocate( bbuffer(tsize*block_buf_nrecs),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'd_p_coupling: failed to allocate bbuffer; error = ',astat
          call endrun
       end if
       allocate( cbuffer(tsize*chunk_buf_nrecs),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'd_p_coupling: failed to allocate cbuffer; error = ',astat
          call endrun
       end if

       if (iam .lt. grid%npes_xy) then
          call block_to_chunk_send_pters(iam+1,blksiz,kmp1,tsize,bpter)
       endif

!$omp parallel do private (j, i, ib, k, m)
!dir$ concurrent
       do j=jfirstxy,jlastxy
!dir$ concurrent
          do i=ifirstxy,ilastxy
             ib = (j-jfirstxy)*(ilastxy-ifirstxy+1) + (i-ifirstxy+1)

             bbuffer(bpter(ib,0)+4:bpter(ib,0)+boff+pcnst) = 0.0_r8

             bbuffer(bpter(ib,0))   = pexy(i,kmp1,j)
             bbuffer(bpter(ib,0)+1) = pelnxy(i,kmp1,j)
             bbuffer(bpter(ib,0)+2) = psxy(i,j)
             bbuffer(bpter(ib,0)+3) = phisxy(i,j)

!dir$ concurrent
             do k=1,km

                bbuffer(bpter(ib,k))   = pexy(i,k,j)
                bbuffer(bpter(ib,k)+1) = pelnxy(i,k,j)
                bbuffer(bpter(ib,k)+2) = u3    (i,k,j)
                bbuffer(bpter(ib,k)+3) = v3    (i,k,j)
                bbuffer(bpter(ib,k)+4) = omgaxy(i,k,j)
                bbuffer(bpter(ib,k)+5) = ptxy(i,j,k) / (D1_0 + zvir*tracer(i,j,k,1))
                bbuffer(bpter(ib,k)+6) = pkxy(i,j,pver+1) / pkzxy(i,j,k) 
                if (do_waccm_phys()) then
                   bbuffer(bpter(ib,k)+7) = frontga(i,k,j)
                   bbuffer(bpter(ib,k)+8) = frontgf(i,k,j)
                   bbuffer(bpter(ib,k)+9) = uzm(k,j)
                end if

                do m=1,pcnst
                   bbuffer(bpter(ib,k)+boff+m) = tracer(i,j,k,m)
                end do

             end do
          end do
       end do

       call t_barrierf('sync_blk_to_chk', grid%commxy)
       call t_startf ('block_to_chunk')
       call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
       call t_stopf  ('block_to_chunk')

!$omp parallel do private (lchnk, ncol, i, k, m, cpter)
chnk_loop2 : &
       do lchnk = begchunk,endchunk
          ncol = phys_state(lchnk)%ncol

          call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

          do i=1,ncol

             phys_state(lchnk)%pint  (i,pver+1) = cbuffer(cpter(i,0))
             phys_state(lchnk)%lnpint(i,pver+1) = cbuffer(cpter(i,0)+1)
             phys_state(lchnk)%ps(i)            = cbuffer(cpter(i,0)+2)
             phys_state(lchnk)%phis(i)          = cbuffer(cpter(i,0)+3)

             do k=1,km

                phys_state(lchnk)%pint  (i,k) = cbuffer(cpter(i,k))
                phys_state(lchnk)%lnpint(i,k) = cbuffer(cpter(i,k)+1)
                phys_state(lchnk)%u     (i,k) = cbuffer(cpter(i,k)+2)
                phys_state(lchnk)%v     (i,k) = cbuffer(cpter(i,k)+3)
                phys_state(lchnk)%omega (i,k) = cbuffer(cpter(i,k)+4)
                phys_state(lchnk)%t     (i,k) = cbuffer(cpter(i,k)+5)
                phys_state(lchnk)%exner (i,k) = cbuffer(cpter(i,k)+6)
                if (do_waccm_phys()) then
                   phys_state(lchnk)%frontga(i,k)  = cbuffer(cpter(i,k)+7)
                   phys_state(lchnk)%frontgf(i,k)  = cbuffer(cpter(i,k)+8)
                   phys_state(lchnk)%uzm(i,k)      = cbuffer(cpter(i,k)+9)
                end if
                
                ! dry type constituents converted from moist to dry at bottom of routine
                do m=1,pcnst
                   phys_state(lchnk)%q(i,k,m) = cbuffer(cpter(i,k)+boff+m)
                end do

             end do
          end do

       end do chnk_loop2

       deallocate(bpter)
       deallocate(bbuffer)
       deallocate(cbuffer)

    endif has_local_map

!------------------------------------------------------
!  Get indices to access O, O2, H, H2, and N species
!------------------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
      call cnst_get_ind('O', ixo)
      call cnst_get_ind('O2', ixo2)
      call cnst_get_ind('H', ixh)
      call cnst_get_ind('H2', ixh2)
      call cnst_get_ind('N', ixn)
    endif
!
! Evaluate derived quantities
!
    call t_startf ('derived_fields')
!$omp parallel do private (lchnk, ncol, i, k, m, qmavl, dqreq, qbot, qbotm1, zvirv, pbuf_chnk, mmrSum_O_O2_H)
    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       do k=1,km
          do i=1,ncol
             phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pint(i,k+1) - phys_state(lchnk)%pint(i,k)
             phys_state(lchnk)%rpdel(i,k) = D1_0/phys_state(lchnk)%pdel(i,k)
             phys_state(lchnk)%pmid (i,k) = D0_5*(phys_state(lchnk)%pint(i,k) + phys_state(lchnk)%pint(i,k+1))
             phys_state(lchnk)%lnpmid(i,k) = log(phys_state(lchnk)%pmid(i,k))
          end do
       end do

! Attempt to remove negative constituents in bottom layer only by moving from next level
! This is a BAB kludge to avoid masses of warning messages for cloud water and ice, since
! the vertical remapping operator currently being used for cam is not strictly monotonic 
! at the endpoints.
       do m=1,pcnst
          do i=1,ncol
             if (phys_state(lchnk)%q(i,pver,m) < qmin(m)) then
! available q in 2nd level
                qmavl = phys_state(lchnk)%q (i,pver-1,m) - qmin(m)
! required q change in bottom level rescaled to mass fraction in 2nd level
                dqreq = (qmin(m) - phys_state(lchnk)%q(i,pver,m))                         &
                      * phys_state(lchnk)%pdel(i,pver) / phys_state(lchnk)%pdel(i,pver-1)
                qbot   = phys_state(lchnk)%q(i,pver  ,m)
                qbotm1 = phys_state(lchnk)%q(i,pver-1,m)
                if (dqreq < qmavl) then
                   phys_state(lchnk)%q(i,pver  ,m) = qmin(m)
                   phys_state(lchnk)%q(i,pver-1,m) = phys_state(lchnk)%q(i,pver-1,m) - dqreq
                   ! Comment out these log messages since they can make the log files so
                   ! large that they're unusable.
                   !if (dqreq>1.e-14_r8 .and. debug_adjust_print) write(iulog,*) 'dpcoup dqreq', m, lchnk, i, qbot, qbotm1, dqreq
                    if (dqreq>qmin(m) .and. dqreq>fraction*qbotm1 .and. debug_adjust_print) &
                                                                write(iulog,*) 'dpcoup dqreq', m, lchnk, i, qbot, qbotm1, dqreq
                else 
                   ! Comment out these log messages since they can make the log files so
                   ! large that they're unusable.
                   !if (debug_adjust_print) write(iulog,*) 'dpcoup cant adjust', m, lchnk, i, qbot, qbotm1, dqreq
                    if (dqreq>qmin(m) .and. debug_adjust_print) write(iulog,*) 'dpcoup cant adjust', m, lchnk, i, &
                         qbot, qbotm1, dqreq
                end if
             end if
          end do
       end do

!-----------------------------------------------------------------------------------------------------------------
! Ensure O2 + O + H (N2) mmr greater than one.  Check for unusually large H2 values and set to lower value
!-----------------------------------------------------------------------------------------------------------------
       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
          do i=1,ncol
             do k=1,pver

                if (phys_state(lchnk)%q(i,k,ixo) < mmrMin) phys_state(lchnk)%q(i,k,ixo) = mmrMin
                if (phys_state(lchnk)%q(i,k,ixo2) < mmrMin) phys_state(lchnk)%q(i,k,ixo2) = mmrMin

                mmrSum_O_O2_H = phys_state(lchnk)%q(i,k,ixo)+phys_state(lchnk)%q(i,k,ixo2)+phys_state(lchnk)%q(i,k,ixh)

                if ((1._r8-mmrMin-mmrSum_O_O2_H) < 0._r8) then
           
                   phys_state(lchnk)%q(i,k,ixo) = phys_state(lchnk)%q(i,k,ixo) * (1._r8 - N2mmrMin) / mmrSum_O_O2_H

                   phys_state(lchnk)%q(i,k,ixo2) = phys_state(lchnk)%q(i,k,ixo2) * (1._r8 - N2mmrMin) / mmrSum_O_O2_H

                   phys_state(lchnk)%q(i,k,ixh) = phys_state(lchnk)%q(i,k,ixh) * (1._r8 - N2mmrMin) / mmrSum_O_O2_H

                endif

                if(phys_state(lchnk)%q(i,k,ixh2) .gt. 6.e-5_r8) then
                   phys_state(lchnk)%q(i,k,ixh2) = 6.e-5_r8
                endif
	     
             end do
          end do
       endif

!-----------------------------------------------------------------------------
! Call physconst_update to compute cpairv, rairv, mbarv, and cappav as constituent dependent variables
! and compute molecular viscosity(kmvis) and conductivity(kmcnd) 
!-----------------------------------------------------------------------------
       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
         call physconst_update(phys_state(lchnk)%q, phys_state(lchnk)%t, &
                                cnst_mw(ixo), cnst_mw(ixo2), cnst_mw(ixh), cnst_mw(ixn), &
                                                           ixo, ixo2, ixh, pcnst, lchnk, ncol)
       endif

!------------------------------------------------------------------------
! Fill local zvirv variable; calculated for WACCM-X
!------------------------------------------------------------------------
       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
         zvirv(:,:) = shr_const_rwv / rairv(:,:,lchnk) -1._r8
       else
         zvirv(:,:) = zvir    
       endif
!
! Compute initial geopotential heights
       call geopotential_t (phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
                            phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
                            phys_state(lchnk)%t     , phys_state(lchnk)%q(1,1,1), rairv(:,:,lchnk), gravit, zvirv, &
                            phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol                )

! Compute initial dry static energy, include surface geopotential
       do k = 1, pver
          do i=1,ncol
             phys_state(lchnk)%s(i,k) = cpairv(i,k,lchnk)*phys_state(lchnk)%t(i,k) &
                                      + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
          end do
       end do

!
! Convert dry type constituents from moist to dry mixing ratio
!
       call set_state_pdry(phys_state(lchnk))    ! First get dry pressure to use for this timestep
       call set_wet_to_dry(phys_state(lchnk))    ! Dynamics had moist, physics wants dry.


! Compute energy and water integrals of input state
       pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
       call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)

    end do
    call t_stopf('derived_fields')

    deallocate (u3)
    deallocate (v3)

!EOC
  end subroutine d_p_coupling
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: p_d_coupling --- convert physics output to dynamics input
!
! !INTERFACE: 
  subroutine p_d_coupling(grid, phys_state, phys_tend, &
                          dyn_in, dtime, zvir, cappa, ptop)

! !USES:
#if ( defined OFFLINE_DYN )
   use metdata,     only: get_met_fields
#endif
 
!-----------------------------------------------------------------------
    implicit none

! Variables ending in xy are xy-decomposition instanciations.

    type(T_FVDYCORE_GRID), intent(in) :: grid ! FV Dynamics grid

! !INPUT PARAMETERS:
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend),  intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(dyn_import_t),  intent(inout)   :: dyn_in

    real(r8), intent(in) :: dtime
    real(r8), intent(in) :: zvir
    real(r8), intent(in) :: cappa
    real(r8), intent(in) :: ptop

! !DESCRIPTION:
!
!   Coupler for converting physics output variables into dynamics input variables
!
! !REVISION HISTORY:
!   00.06.01   Boville    Creation
!   01.06.08   AAM        Compactified
!   01.07.13   AAM        Some support for multi-2D decompositions
!   02.03.01   Worley     Support for nontrivial physics remapping
!   02.08.06   Sawyer     T3 added -- updated to current temperature
!   05.07.12   Sawyer     Added dyn_state as argument
!   05.09.23   Sawyer     Transitioned to XY decomposition vars. only
!   05.10.31   Sawyer     Replaced dyn_state with dyn_interface
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

! Variables from the dynamics import container

    real(r8), pointer :: psxy(:,:)
    real(r8), pointer :: u3sxy(:,:,:)
    real(r8), pointer :: v3sxy(:,:,:)
    real(r8), pointer :: t3xy(:,:,:)                  !  Temperature
    real(r8), pointer :: ptxy(:,:,:)                  !  Virt. pot. temp.
    real(r8), pointer :: tracer(:,:,:,:)              !  Constituents

    real(r8), pointer :: pexy(:,:,:)
    real(r8), pointer :: delpxy(:,:,:)
    real(r8), pointer :: pkxy(:,:,:)
    real(r8), pointer :: pkzxy(:,:,:)

! Local workspace

    real(r8):: dudtxy(grid%ifirstxy:grid%ilastxy,&
                      grid%km,grid%jfirstxy:grid%jlastxy)
    real(r8):: dvdtxy(grid%ifirstxy:grid%ilastxy,&
                      grid%km,grid%jfirstxy:grid%jlastxy)
    real(r8):: dummy_pelnxy(grid%ifirstxy:grid%ilastxy,grid%km+1, &
                            grid%jfirstxy:grid%jlastxy)

    integer :: i, ib, k, m, j, lchnk  ! indices
    integer :: ncol                   ! number of columns in current chunk
    integer :: lats(pcols)            ! array of latitude indices
    integer :: lons(pcols)            ! array of longitude indices
    integer :: blksiz                 ! number of columns in 2D block
    integer :: tsize                  ! amount of data per grid point passed to physics
    integer, allocatable, dimension(:,:) :: bpter
                                     ! offsets into block buffer for unpacking data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for packing data
    integer :: iqa, iqb, iqc, iqd, mq     ! used for tracer transpose grouping

    real(r8) :: dt5
    real(r8), allocatable, dimension(:) :: &
       bbuffer, cbuffer               ! transpose buffers
#if (! defined SPMD)
    integer  :: block_buf_nrecs = 0
    integer  :: chunk_buf_nrecs = 0
    logical  :: local_dp_map=.true. 
#endif
    integer  :: im, jm, km, ng_d, ng_s, iam
    integer  :: ifirstxy, ilastxy, jfirstxy, jlastxy 
    integer  :: jfirst, jlast, kfirst, klast

! Pull the variables out of the dynamics export container

    psxy    => dyn_in%ps
    u3sxy   => dyn_in%u3s
    v3sxy   => dyn_in%v3s
    t3xy    => dyn_in%t3
    ptxy    => dyn_in%pt
    tracer  => dyn_in%tracer

    pexy    => dyn_in%pe
    delpxy  => dyn_in%delp
    pkxy    => dyn_in%pk
    pkzxy   => dyn_in%pkz    

    im   = grid%im
    jm   = grid%jm
    km   = grid%km

    ifirstxy = grid%ifirstxy
    ilastxy  = grid%ilastxy
    jfirstxy = grid%jfirstxy
    jlastxy  = grid%jlastxy

    jfirst   = grid%jfirst
    jlast    = grid%jlast
    kfirst   = grid%kfirst
    klast    = grid%klast

    ng_d     = grid%ng_d
    ng_s     = grid%ng_s

    iam      = grid%iam

!---------------------------End Local workspace-------------------------

#if ( defined OFFLINE_DYN )
!
! set the dyn flds to offline meteorological data
!
      call get_met_fields( phys_state, phys_tend, dtime )
#endif
! -------------------------------------------------------------------------
! Copy temperature, tendencies and constituents to dynamics data structures
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
! Copy onto xy decomposition, then transpose to yz decomposition
! -------------------------------------------------------------------------

       if (local_dp_map) then

!$omp parallel do private(lchnk, i, k, ncol, m, lons, lats)

          do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)
             call get_lon_all_p(lchnk, ncol, lons)
             call get_lat_all_p(lchnk, ncol, lats)

             do k = 1, km
                do i = 1, ncol
                   dvdtxy(lons(i),k,lats(i)) = phys_tend(lchnk)%dvdt(i,k)
                   dudtxy(lons(i),k,lats(i)) = phys_tend(lchnk)%dudt(i,k)
                   ptxy  (lons(i),lats(i),k) = phys_state(lchnk)%t(i,k)
                   delpxy(lons(i),lats(i),k) = phys_state(lchnk)%pdel(i,k)
                enddo
             enddo

             do m=1,pcnst
                do k=1,km
                   do i=1,ncol
                      tracer(lons(i),lats(i),k,m) = &
                         phys_state(lchnk)%q(i,k,m)
                   end do
                end do
             end do

          enddo

       else

          tsize = 4 + pcnst

          blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
          allocate(bpter(blksiz,0:km))
          allocate(bbuffer(tsize*block_buf_nrecs))
          allocate(cbuffer(tsize*chunk_buf_nrecs))

!$omp parallel do private (lchnk, ncol, i, k, m, cpter)
          do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)

             call chunk_to_block_send_pters(lchnk,pcols,km+1,tsize,cpter)

             do i=1,ncol
                cbuffer(cpter(i,0):cpter(i,0)+3+pcnst) = 0.0_r8
             end do

!dir$ concurrent
             do k=1,km
!dir$ concurrent
                do i=1,ncol

                   cbuffer(cpter(i,k))   = phys_tend(lchnk)%dvdt(i,k)
                   cbuffer(cpter(i,k)+1) = phys_tend(lchnk)%dudt(i,k)
                   cbuffer(cpter(i,k)+2) = phys_state(lchnk)%t(i,k)
                   cbuffer(cpter(i,k)+3) = phys_state(lchnk)%pdel(i,k)

                   do m=1,pcnst
                      cbuffer(cpter(i,k)+3+m) = phys_state(lchnk)%q(i,k,m)
                   end do

                end do
  
             end do

          end do

          call t_barrierf('sync_chk_to_blk', grid%commxy)
          call t_startf ('chunk_to_block')
          call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
          call t_stopf  ('chunk_to_block')

          if (iam .lt. grid%npes_xy) then
             call chunk_to_block_recv_pters(iam+1,blksiz,km+1,tsize,bpter)
          endif

!$omp parallel do private (j, i, ib, k, m)
!dir$ concurrent
          do j=jfirstxy,jlastxy
!dir$ concurrent
             do k=1,km
!dir$ concurrent
                do i=ifirstxy,ilastxy
                   ib = (j-jfirstxy)*(ilastxy-ifirstxy+1) + (i-ifirstxy+1)

                   dvdtxy(i,k,j) = bbuffer(bpter(ib,k))
                   dudtxy(i,k,j) = bbuffer(bpter(ib,k)+1)
                   ptxy  (i,j,k) = bbuffer(bpter(ib,k)+2)
                   delpxy(i,j,k) = bbuffer(bpter(ib,k)+3)

                   do m=1,pcnst
                      tracer(i,j,k,m) = bbuffer(bpter(ib,k)+3+m)
                   end do

                enddo
             enddo
          enddo

          deallocate(bpter)
          deallocate(bbuffer)
          deallocate(cbuffer)

       endif

! WS: 02.08.06: Update t3 to temperature
!$omp parallel do private(i,j,k)
!dir$ concurrent
       do k=1,km
          do j = jfirstxy,jlastxy
             do i = ifirstxy,ilastxy
                t3xy(i,j,k) = ptxy(i,j,k)
             enddo
          enddo
       enddo

! -------------------------------------------------------------------------
! Update u3s and v3s from tendencies dudt and dvdt.
! -------------------------------------------------------------------------
       dt5 = D0_5*dtime

       call t_barrierf('sync_uv3s_update', grid%commxy)
       call t_startf('uv3s_update')
       if (iam .lt. grid%npes_xy) then
          call uv3s_update( grid, dudtxy, u3sxy, dvdtxy, v3sxy, dt5 )
       end if  ! (iam .lt. grid%npes_xy)
       call t_stopf('uv3s_update')

! -------------------------------------------------------------------------
! Compute pt, q3, pe, delp, ps, peln, pkz and pk.
! For 2-D decomposition, delp is transposed to delpxy, pexy is computed
!  from delpxy (and ptop), and pexy is transposed back to pe.
! Note that pt, q3, delp and pe are input parameters as well.
! -------------------------------------------------------------------------
    call t_barrierf('sync_p_d_adjust', grid%commxy)
    call t_startf ('p_d_adjust')
    if (iam .lt. grid%npes_xy) then
       call p_d_adjust(grid, tracer, dummy_pelnxy, pkxy, pkzxy, zvir,  cappa, &
                       delpxy, ptxy, pexy, psxy, ptop)
    end if  ! (iam .lt. grid%npes_xy)
    call t_stopf  ('p_d_adjust')

!EOC
  end subroutine p_d_coupling
!-----------------------------------------------------------------------
end module dp_coupling
