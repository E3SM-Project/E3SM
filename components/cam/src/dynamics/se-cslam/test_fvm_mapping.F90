module test_fvm_mapping
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use fvm_control_volume_mod, only: fvm_struct
  use cam_history,            only: outfld
  use physconst,              only: pi
  use dimensions_mod,         only: qsize_condensate_loading,qsize_condensate_loading_idx
  use dimensions_mod,         only: np, nelemd, nlev, npsq, ntrac
  use element_mod,            only: element_t
  implicit none
  private

  real(r8), parameter, private :: deg2rad = pi/180.0_r8
  real(r8), parameter, private :: psurf_moist = 100000.0_r8 !moist surface pressure
  integer,  parameter, private :: num_tracer = 13
  integer,  parameter, private :: cl_idx  = 12
  integer,  parameter, private :: cl2_idx = 13
  real(r8), parameter, private :: cly_constant = 4.e-6_r8

  public :: test_mapping_overwrite_tendencies, test_mapping_addfld
  public :: test_mapping_output_mapped_tendencies, test_mapping_overwrite_dyn_state
  public :: test_mapping_output_phys_state
contains

  subroutine test_mapping_addfld
#ifdef debug_coupling
    use cam_history,        only: addfld, add_default, horiz_only, register_vector_field
    use constituents,       only: cnst_get_ind,cnst_name
    character(LEN=128) :: name
    integer :: nq,m_cnst

    name = 'd2p_u_gll'
    call addfld(trim(name),   (/ 'lev' /),  'I','m/2','Exact zonal wind on GLL grid',gridname='GLL')
    call add_default (trim(name), 1, ' ')

    name = 'd2p_v_gll'
    call addfld(trim(name),   (/ 'lev' /),  'I','m/2','Exact meridional wind on GLL grid',gridname='GLL')
    call add_default (trim(name), 1, ' ')

    name = 'd2p_scalar_gll'
    call addfld(trim(name),   (/ 'lev' /),  'I','','Exact scalar on GLL grid',gridname='GLL')
    call add_default (trim(name), 1, ' ')

    name = 'd2p_u'
    call addfld(trim(name),   (/ 'lev' /),  'I','m/2','Zonal wind mapped to physics grid')
    call add_default (trim(name), 1, ' ')

    name = 'd2p_u_err'
    call addfld(trim(name),   (/ 'lev' /),  'I','m/2','Error in zonal wind mapped to physics grid')
    call add_default (trim(name), 1, ' ')

    name = 'd2p_v_err'
    call addfld(trim(name),   (/ 'lev' /),  'I','m/2','Error in meridional wind mapped to physics grid')
    call add_default (trim(name), 1, ' ')

    name = 'd2p_v'
    call addfld(trim(name),   (/ 'lev' /),  'I','m/s','Meridional wind mapped to physics grid')
    call add_default (trim(name), 1, ' ')

    name = 'd2p_scalar'
    call addfld(trim(name),   (/ 'lev' /),  'I','','Scalar mapped to physics grid')
    call add_default (trim(name), 1, ' ')

    name = 'd2p_scalar_err'
    call addfld(trim(name),   (/ 'lev' /),  'I','','Error in scalar mapped to physics grid')
    call add_default (trim(name), 1, ' ')

    do nq=2,qsize_condensate_loading
      m_cnst = qsize_condensate_loading_idx(nq)
      name = 'f2p_'//trim(cnst_name(m_cnst))//'_fvm'
      call addfld(trim(name),   (/ 'lev' /),  'I','','Exact water tracer on fvm grid',gridname='FVM')
      call add_default (trim(name), 1, ' ')
      name = 'f2p_'//trim(cnst_name(m_cnst))//'_err'
      call addfld(trim(name),   (/ 'lev' /),  'I','','Error in water tracer on physics grid (mapped from fvm grid)')
      call add_default (trim(name), 1, ' ')
      name = 'f2p_'//trim(cnst_name(m_cnst))//''
      call addfld(trim(name),   (/ 'lev' /),  'I','','Water tracer on physics grid (mapped from fvm grid')
      call add_default (trim(name), 1, ' ')
      !
      ! physgrid to gll (condensate loading tracers)
      !
      name = 'p2d_'//trim(cnst_name(m_cnst))//''
      call addfld(trim(name),   (/ 'lev' /),  'I','','Water tracer on physics grid')
      call add_default (trim(name), 1, ' ')
      name = 'p2d_'//trim(cnst_name(m_cnst))//'_gll'
      call addfld(trim(name),   (/ 'lev' /),  'I','','Water tracer on GLL grid',gridname='GLL')
      call add_default (trim(name), 1, ' ')
      name = 'p2d_'//trim(cnst_name(m_cnst))//'_err_gll'
      call addfld(trim(name),   (/ 'lev' /),  'I','','Error in water tracer mapped to GLL grid',gridname='GLL')
      call add_default (trim(name), 1, ' ')
      !
      ! physgrid to fvm (condensate loading tracers)
      !
      name = 'p2f_'//trim(cnst_name(m_cnst))//''
      call addfld(trim(name),   (/ 'lev' /),  'I','','Water tracer on physics grid')
      call add_default (trim(name), 1, ' ')
      name = 'p2f_'//trim(cnst_name(m_cnst))//'_fvm'
      call addfld(trim(name),   (/ 'lev' /),  'I','','Water tracer on FVM grid',gridname='FVM')
      call add_default (trim(name), 1, ' ')
      name = 'p2f_'//trim(cnst_name(m_cnst))//'_err_fvm'
      call addfld(trim(name),   (/ 'lev' /),  'I','','Error in water tracer mapped to FVM grid',gridname='FVM')
      call add_default (trim(name), 1, ' ')
    end do
    !
    ! temperature tendency
    !
    name = 'p2d_ptend'
    call addfld(trim(name),   (/ 'lev' /),  'I','','T tendency on physics grid')
    call add_default (trim(name), 1, ' ')
    name = 'p2d_ptend_gll'
    call addfld(trim(name),   (/ 'lev' /),  'I','','T tendency on GLL grid',gridname='GLL')
    call add_default (trim(name), 1, ' ')
    name = 'p2d_ptend_err_gll'
    call addfld(trim(name),   (/ 'lev' /),  'I','','Error in T tendency mapped to GLL grid',gridname='GLL')
    call add_default (trim(name), 1, ' ')

    call addfld('p2d_u',   (/ 'lev' /),  'I','m/2','Zonal wind on physics grid')
    call add_default ('p2d_u', 1, ' ')
    call addfld('p2d_v',   (/ 'lev' /),  'I','m/2','Meridional wind on physics grid')
    call add_default ('p2d_v', 1, ' ')
    call addfld('p2d_u_gll',   (/ 'lev' /),  'I','m/2','Zonal wind on physics grid',gridname='GLL')
    call add_default ('p2d_u_gll', 1, ' ')
    call addfld('p2d_v_gll',   (/ 'lev' /),  'I','m/2','Meridional wind on physics grid',gridname='GLL')
    call add_default ('p2d_v_gll', 1, ' ')
    call addfld('p2d_u_gll_err',   (/ 'lev' /),  'I','m/2','Error in zonal wind interpolation to GLL grid',gridname='GLL')
    call add_default ('p2d_u_gll_err', 1, ' ')
    call addfld('p2d_v_gll_err',   (/ 'lev' /),  'I','m/2','Error in meridional wind interpolation to GLL grid',&
         gridname='GLL')
    call add_default ('p2d_v_gll_err', 1, ' ')

!      name = 'phys2dyn_'//trim(cnst_name(m_cnst))//'_physgrid'
!      call outfld(trim(name),phys_state%q(:ncols,:,m_cnst),ncols,lchnk)
#endif
  end subroutine test_mapping_addfld

  subroutine test_mapping_overwrite_tendencies(phys_state,phys_tend,ncols,lchnk,q_prev)
    use dimensions_mod,         only: fv_nphys
    use constituents,           only: cnst_get_ind,pcnst,cnst_name
    use physics_types,  only: physics_state, physics_tend
    type(physics_state), intent(inout) :: phys_state
    type(physics_tend),  intent(inout) :: phys_tend
    real(r8), dimension(:,:,:), intent(inout) :: q_prev
!    type(fvm_struct), intent(inout):: fvm(:)
    integer,          intent(in)   :: ncols,lchnk
#ifdef debug_coupling
    integer :: icol,k
    character(LEN=128) :: name
    integer :: m_cnst, nq

    q_prev = 0.0_r8
    phys_state%pdel(1:ncols,:) = phys_state%pdeldry(1:ncols,:) !make sure there is no conversion from wet to dry
    do nq=2,qsize_condensate_loading
      m_cnst = qsize_condensate_loading_idx(nq)
      do icol=1,ncols
        do k=1,num_tracer
          phys_state%q(icol,k,m_cnst)   = test_func(phys_state%lat(icol), phys_state%lon(icol), k, k)
        end do
      enddo
      name = 'p2f_'//trim(cnst_name(m_cnst))//''
      call outfld(trim(name),phys_state%q(:ncols,:,m_cnst),ncols,lchnk)
      name = 'p2d_'//trim(cnst_name(m_cnst))//''
      call outfld(trim(name),phys_state%q(:ncols,:,m_cnst),ncols,lchnk)
    end do

    do icol=1,ncols
      do k=1,num_tracer
        phys_tend%dudt(icol,k) = test_func(phys_state%lat(icol), phys_state%lon(icol), k, k)
        phys_tend%dvdt(icol,k) = test_func(phys_state%lat(icol), phys_state%lon(icol), k, k)
        phys_tend%dtdt(icol,k) = test_func(phys_state%lat(icol), phys_state%lon(icol), k, k)
      end do
    enddo
    name = 'p2d_u'
    call outfld(trim(name),phys_tend%dudt(:ncols,:),ncols,lchnk)
    name = 'p2d_v'
    call outfld(trim(name),phys_tend%dvdt(:ncols,:),ncols,lchnk)
    name = 'p2d_ptend'
    call outfld(trim(name),phys_tend%dtdt(:ncols,:),ncols,lchnk)


!    do icol=1,ncols
!      do k=1,nlev
!        phys_tend%dudt(icol,k) = 0.0_r8
!        phys_tend%dvdt(icol,k) = 0.0_r8
!      end do
!    enddo
#endif
  end subroutine test_mapping_overwrite_tendencies

  subroutine test_mapping_output_mapped_tendencies(fvm,elem,nets,nete,tl_f,tl_qdp)
    use dimensions_mod,         only: fv_nphys,nlev,nc
    use constituents,           only: cnst_get_ind,cnst_name
    integer,          intent(in)   :: nets,nete,tl_f,tl_qdp
    type(fvm_struct), intent(inout):: fvm(nets:nete)
    type(element_t),  intent(inout):: elem(nets:nete)             ! pointer to dyn_out element array
#ifdef debug_coupling
    integer :: ie,i,j,k
    character(LEN=128) :: name
    integer :: nq,m_cnst
    real(r8) :: diff(nc,nc,nlev)

    do ie = nets,nete
      call outfld('p2d_u_gll', RESHAPE(elem(ie)%derived%fm(:,:,1,:),(/npsq,nlev/)), npsq, ie)
      call outfld('p2d_v_gll', RESHAPE(elem(ie)%derived%fm(:,:,2,:),(/npsq,nlev/)), npsq, ie)
      call outfld('p2d_ptend_gll', RESHAPE(elem(ie)%derived%ft(:,:,:),(/npsq,nlev/)), npsq, ie)
      do k=1,num_tracer
        do j=1,np
          do i=1,np
            elem(ie)%derived%fm(i,j,1,k)   = elem(ie)%derived%fm(i,j,1,k) -&
                 test_func(elem(ie)%spherep(i,j)%lat, elem(ie)%spherep(i,j)%lon, k,k)
            elem(ie)%derived%fm(i,j,2,k)   = elem(ie)%derived%fm(i,j,2,k) - &
                 test_func(elem(ie)%spherep(i,j)%lat, elem(ie)%spherep(i,j)%lon, k,k)
            elem(ie)%derived%ft(i,j,k)   = elem(ie)%derived%ft(i,j,k) - &
                 test_func(elem(ie)%spherep(i,j)%lat, elem(ie)%spherep(i,j)%lon, k,k)
          end do
        end do
      end do
      call outfld('p2d_u_gll_err'    , RESHAPE(elem(ie)%derived%fm(:,:,1,:),(/npsq,nlev/)), npsq, ie)
      call outfld('p2d_v_gll_err'    , RESHAPE(elem(ie)%derived%fm(:,:,2,:),(/npsq,nlev/)), npsq, ie)
      call outfld('p2d_ptend_err_gll', RESHAPE(elem(ie)%derived%ft(:,:,:),(/npsq,nlev/)), npsq, ie)
      elem(ie)%derived%ft(:,:,:)   = 0.0_r8
    end do

    do ie = nets,nete
      do nq=2,qsize_condensate_loading
        m_cnst = qsize_condensate_loading_idx(nq)
        name = 'p2d_'//trim(cnst_name(m_cnst))//'_gll'
        call outfld(TRIM(name), RESHAPE(elem(ie)%derived%fq(:,:,:,nq),(/npsq,nlev/)), npsq, ie)
        !        call outfld(trim(name),&
        !             RESHAPE(fvm(ie)%fc(1:nc,1:nc,:,m_cnst),&
        !             (/nc*nc,nlev/)),nc*nc,ie)
        do k=1,num_tracer
          do j=1,np
            do i=1,np
              elem(ie)%derived%fq(i,j,k,nq) = elem(ie)%derived%fq(i,j,k,nq)-&
                   test_func(elem(ie)%spherep(i,j)%lat, elem(ie)%spherep(i,j)%lon, k, k)
            end do
          end do
        end do
        name = 'p2d_'//trim(cnst_name(m_cnst))//'_err_gll'
        call outfld(TRIM(name), RESHAPE(elem(ie)%derived%fq(:,:,:,nq),(/npsq,nlev/)), npsq, ie)
      end do
      if (ntrac>0) then
        do nq=2,qsize_condensate_loading
          m_cnst = qsize_condensate_loading_idx(nq)
          name = 'p2f_'//trim(cnst_name(m_cnst))//'_fvm'
          !
          ! cly
          !
!          k=num_tracer+1
!          fvm(ie)%fc(1:nc,1:nc,k,:) = fvm(ie)%fc(1:nc,1:nc,cl_idx,:)+&
!                                    2.0_r8*fvm(ie)%fc(1:nc,1:nc,cl2_idx,:)
          call outfld(trim(name),&
               RESHAPE(fvm(ie)%fc(1:nc,1:nc,:,m_cnst),&
               (/nc*nc,nlev/)),nc*nc,ie)
          do k=1,num_tracer
            do j=1,nc
              do i=1,nc
                fvm(ie)%fc(i,j,k,m_cnst) = fvm(ie)%fc(i,j,k,m_cnst)-&
                     test_func(fvm(ie)%center_cart(i,j)%lat,fvm(ie)%center_cart(i,j)%lon, k, k)
              end do
            end do
          end do
          name = 'p2f_'//trim(cnst_name(m_cnst))//'_err_fvm'
          call outfld(TRIM(name), RESHAPE(fvm(ie)%fc(:,:,:,m_cnst),(/nc*nc,nlev/)), nc*nc, ie)

        end do
      endif
    end do
#endif
  end subroutine test_mapping_output_mapped_tendencies

  subroutine test_mapping_overwrite_dyn_state(elem,fvm)
    use fvm_control_volume_mod, only: fvm_struct,n0_fvm
    use constituents,           only: cnst_name
    use dimensions_mod,         only: nc,nhc
    use hybrid_mod,             only: get_loop_ranges, hybrid_t,config_thread_region
!    use fvm_mod,                only: fill_halo_fvm_noprealloc
    use parallel_mod,           only: par
    type (fvm_struct), intent(inout)    :: fvm(:)
    type(element_t), intent(inout) :: elem(:)             ! pointer to dyn_out element array
#ifdef debug_coupling
    integer            :: i,j,k,ie,nq,m_cnst
    character(LEN=128) :: name
    integer            :: nets,nete
    type(hybrid_t)                                                   :: hybrid

    hybrid = config_thread_region(par,'serial')
    call get_loop_ranges(hybrid,ibeg=nets,iend=nete)
    do ie=nets,nete
      do nq=2,qsize_condensate_loading        
        m_cnst = qsize_condensate_loading_idx(nq)
        name = 'f2p_'//trim(cnst_name(m_cnst))//'_fvm'
        do k=1,num_tracer
          do j=1,nc
            do i=1,nc
              fvm(ie)%c(i,j,k,m_cnst,:) = test_func(fvm(ie)%center_cart(i,j)%lat,fvm(ie)%center_cart(i,j)%lon, k, k)
            end do
          end do
        end do
        !
        ! cly
        !
        k=num_tracer+1
        do j=1,nc
          do i=1,nc
            fvm(ie)%c(i,j,k,m_cnst,:) = test_func(fvm(ie)%center_cart(i,j)%lat,fvm(ie)%center_cart(i,j)%lon, k,cl_idx)+&
                                 2.0_r8*test_func(fvm(ie)%center_cart(i,j)%lat,fvm(ie)%center_cart(i,j)%lon, k,cl2_idx)
          end do
        end do
        call outfld(TRIM(name), RESHAPE(fvm(ie)%c(1:nc,1:nc,:,m_cnst,1),(/nc*nc,nlev/)), nc*nc, ie)
      end do

      elem(ie)%state%Qdp(:,:,:,:,:)   = 0.0_r8 !for testing the p2d map
      do k=1,num_tracer
        do j=1,np
          do i=1,np
            elem(ie)%state%v(i,j,1,k,:)   = test_func(elem(ie)%spherep(i,j)%lat, elem(ie)%spherep(i,j)%lon, k, k )
            elem(ie)%state%v(i,j,2,k,:)   = test_func(elem(ie)%spherep(i,j)%lat, elem(ie)%spherep(i,j)%lon, k, k)
          end do
        end do
      end do
      do k=1,num_tracer
        do j=1,np
          do i=1,np
            elem(ie)%derived%omega(i,j,k) = test_func(elem(ie)%spherep(i,j)%lat, elem(ie)%spherep(i,j)%lon, k, k)
          end do
        end do
      end do
      call outfld('d2p_scalar_gll', RESHAPE(elem(ie)%derived%omega(:,:,:)  ,(/npsq,nlev/)), npsq, ie)
      call outfld('d2p_u_gll', RESHAPE(elem(ie)%state%v(:,:,1,:,1),(/npsq,nlev/)), npsq, ie)
      call outfld('d2p_v_gll', RESHAPE(elem(ie)%state%v(:,:,2,:,1),(/npsq,nlev/)), npsq, ie)
    end do
    !
    ! do boundary exchange (this call should be indentical to call in prim_driver)
    !
!    call fill_halo_fvm_noprealloc(elem,fvm,hybrid,nets,nete,n0_fvm,nhc,1,nlev)!xxx nhr chould be a function of interp_method
#endif
  end subroutine test_mapping_overwrite_dyn_state

  subroutine test_mapping_output_phys_state(phys_state,fvm)
    use physics_types, only: physics_state
    use ppgrid,        only: begchunk, endchunk, pver, pcols
    use constituents,  only: cnst_get_ind,cnst_name
    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(fvm_struct), pointer:: fvm(:)
#ifdef debug_coupling
    integer            :: lchnk, ncol,k,icol,m_cnst,nq,ie
    character(LEN=128) :: name

    do ie=1,nelemd
      fvm(ie)%c(:,:,:,:,:) = 0.0_r8
    end do

    do lchnk = begchunk, endchunk
      call outfld('d2p_scalar', phys_state(lchnk)%omega(1:pcols,1:pver), pcols, lchnk)
      call outfld('d2p_u', phys_state(lchnk)%U(1:pcols,1:pver), pcols, lchnk)
      call outfld('d2p_v', phys_state(lchnk)%V(1:pcols,1:pver), pcols, lchnk)
      if (ntrac>0) then
        do nq=2,qsize_condensate_loading
          m_cnst = qsize_condensate_loading_idx(nq)
          name = 'f2p_'//trim(cnst_name(m_cnst))
          !
          ! cly
          !
          !phys_state(lchnk)%q(1:pcols,num_tracer+1,m_cnst)=phys_state(lchnk)%q(1:pcols,cl_idx,m_cnst)+&
          !      2.0_r8*phys_state(lchnk)%q(1:pcols,12,m_cnst)
          call outfld(TRIM(name), phys_state(lchnk)%q(1:pcols,1:pver,m_cnst), pcols, lchnk)
          k=num_tracer+1
          do icol=1,phys_state(lchnk)%ncol
            phys_state(lchnk)%q(icol,k,m_cnst) = phys_state(lchnk)%q(icol,cl_idx,m_cnst)+&
                                          2.0_r8*phys_state(lchnk)%q(icol,cl2_idx,m_cnst)-&
                                                 cly_constant
          end do
          do k=1,num_tracer
            do icol=1,phys_state(lchnk)%ncol
              phys_state(lchnk)%q(icol,k,m_cnst) = phys_state(lchnk)%q(icol,k,m_cnst)&
                   -test_func(phys_state(lchnk)%lat(icol), phys_state(lchnk)%lon(icol), k,k)
            end do
          enddo
          name = 'f2p_'//trim(cnst_name(m_cnst))//'_err'
          call outfld(TRIM(name), phys_state(lchnk)%q(1:pcols,1:pver,m_cnst), pcols, lchnk)
          phys_state(lchnk)%q(1:pcols,1:pver,m_cnst) = 0.0_r8
        end do
      end if
    end do


    do lchnk = begchunk, endchunk
      do k=1,nlev
        do icol=1,phys_state(lchnk)%ncol
          phys_state(lchnk)%U(icol,k) = phys_state(lchnk)%U(icol,k)&
               -test_func(phys_state(lchnk)%lat(icol), phys_state(lchnk)%lon(icol), k, 9)
          phys_state(lchnk)%V(icol,k) = phys_state(lchnk)%V(icol,k)&
               -test_func(phys_state(lchnk)%lat(icol), phys_state(lchnk)%lon(icol), k,10)
        end do
      enddo
      name = 'd2p_u_err'
      call outfld(trim(name),phys_state(lchnk)%U(:pcols,:),pcols,lchnk)
      name = 'd2p_v_err'
      call outfld(trim(name),phys_state(lchnk)%V(:pcols,:),pcols,lchnk)
      do k=1,num_tracer
        do icol=1,phys_state(lchnk)%ncol
          phys_state(lchnk)%omega(icol,k) = phys_state(lchnk)%omega(icol,k)&
               -test_func(phys_state(lchnk)%lat(icol), phys_state(lchnk)%lon(icol), k,k)
        end do
      end do
      name = 'd2p_scalar_err'
      call outfld(trim(name),phys_state(lchnk)%omega(:pcols,:),pcols,lchnk)
    end do
#endif
  end subroutine test_mapping_output_phys_state

!  subroutine test_mapping_overwrite_state(phys_tend,nets,nete)
!#ifdef debug_coupling
!    phys_tend(lchnk)%dtdt(icol,ilyr)   = 0!test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,9)
!    phys_tend(lchnk)%dudt(icol,ilyr)   = 0!test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,12)
!    phys_tend(lchnk)%dvdt(icol,ilyr)   = 0!test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,13)
!    q_prev(icol, ilyr, 2:pcnst, lchnk) = 0.0D0
!    do m=2,pcnst
!      phys_state(lchnk)%Q(icol,ilyr,m)=0!test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,m)
!    end do
!#endif
!  end subroutine test_mapping_overwrite_state

#ifdef debug_coupling
  function test_func(lat_in, lon_in, k, funcnum) result(fout)
    use hycoef,        only: hyai, hybi, hyam, hybm, ps0
    use shr_sys_mod,   only: shr_sys_flush
    use cam_abortutils, only: endrun
    real(r8), intent(in) :: lon_in
    real(r8), intent(in) :: lat_in
    integer,  intent(in) :: k
    integer,  intent(in) :: funcnum
    real(r8)             :: fout
    real(r8)             :: lon1,lat1,R0,Rg1,Rg2,lon2,lat2,cl,cl2
    real(r8)             :: eta_c

    real(r8)             :: radius      = 10.0_r8 ! radius of the perturbation
    real(r8)             :: perturb_lon = 20.0_r8 ! longitudinal position, 20E
    real(r8)             :: perturb_lat = 40.0_r8 ! latitudinal position, 40N
    real(r8)             :: cos_tmp, sin_tmp, eta
    real(r8)             :: u_wind, v_wind, lat, lon, u_tmp, v_tmp
    real(r8)             :: rotation_angle
    real(r8)             :: det,r,k1,k2
    real(r8), parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164_r8
    real(r8), parameter :: half_pi = pi*0.5_r8
    real(r8), parameter :: degrees_to_radians = pi/180.0_r8
    real(r8), parameter :: k1_lat_center =   20.d0*degrees_to_radians
    real(r8), parameter :: k1_lon_center =  300.d0*degrees_to_radians

    lon = lon_in
    lat = lat_in


    select case(funcnum)
    case(1)
      !
      !   Non-smooth scalar field (slotted cylinder)
      !
      R0 = 0.5_r8
      lon1 = 4.0_r8 * PI / 5.0_r8
      lat1 = 0.0_r8
      Rg1 = acos(sin(lat1)*sin(lat)+cos(lat1)*cos(lat)*cos(lon-lon1))
      lon2 = 6.0_r8 * PI / 5.0_r8
      lat2 = 0.0_r8
      Rg2 = acos(sin(lat2)*sin(lat)+cos(lat2)*cos(lat)*cos(lon-lon2))

      if ((Rg1 <= R0) .AND. (abs(lon-lon1) >= R0/6)) then
        fout = 2.0_r8
      elseif ((Rg2 <= R0) .AND. (abs(lon-lon2) >= R0/6)) then
        fout = 2.0_r8
      elseif ((Rg1 <= R0) .AND. (abs(lon-lon1) < R0/6) &
           .AND. (lat-lat1 < -5.0_r8*R0/12.0_r8)) then
        fout = 2.0_r8
      elseif ((Rg2 <= R0) .AND. (abs(lon-lon2) < R0/6) &
           .AND. (lat-lat2 > 5.0_r8*R0/12.0_r8)) then
        fout = 2.0_r8
      else
        fout = 1.0_r8
      endif
    case(2)
      !
      ! Smooth Gaussian "ball"
      !
      R0    = 10.0_r8           ! radius of the perturbation
      lon1  = 20.0_r8*deg2rad   ! longitudinal position, 20E
      lat1  = 40.0_r8 *deg2rad  ! latitudinal position, 40N
      eta_c = 0.6_r8
      sin_tmp = SIN(lat1)*SIN(lat)
      cos_tmp = COS(lat1)*COS(lat)
      Rg1 = ACOS( sin_tmp + cos_tmp*COS(lon-lon1) )    ! great circle distance
      eta =  (hyam(k)*ps0 + hybm(k)*psurf_moist)/psurf_moist
      fout = EXP(- ((Rg1*R0)**2 + ((eta-eta_c)/0.1_r8)**2))
      if (ABS(fout) < 1.0E-8_r8) then
        fout = 0.0_r8
      end IF
    case(3)
      !
      !
      !
      fout = 0.5_r8 * ( tanh( 3.0_r8*abs(lat)-pi ) + 1.0_r8)
    case(4)
      fout = 2.0_r8+cos(5.0_r8+40*lon)!1.0e-8_r8
    case(5)
      !
      ! approximately Y^2_2 spherical harmonic
      !
      fout = sin(lon)*cos(40*lat)!1.0e-8_r8
!xxx      fout = 0.5_r8 + 0.5_r8*(cos(lat)*cos(lat)*cos(2.0_r8*lon))
    case(6)
      !
      ! approximately Y32_16 spherical harmonic
      !
      fout = 0.5_r8 + 0.5_r8*(cos(16*lon)*(sin(2_r8*lat)**16))
    case(7)
      fout = 2.0_r8 + lat
    case(8)
      fout = 2.0_r8 + cos(lon)
    case(9)
      rotation_angle = 45.0_r8*pi/180.0_r8
      CALL regrot(lon_in,lat_in,lon,lat,0.0_r8,-0.5_r8*pi+rotation_angle,1)
      call Rossby_Haurwitz (lon, lat,u_wind, v_wind)
      CALL turnwi(u_wind,v_wind,u_tmp,v_tmp,lon_in,lat_in,lon,lat,0.0_r8,-0.5_r8*pi+rotation_angle,-1)
      fout = u_tmp
    case(10)
      rotation_angle = 45.0_r8*pi/180.0_r8
      CALL regrot(lon_in,lat_in,lon,lat,0.0_r8,-0.5_r8*pi+rotation_angle,1)
      call Rossby_Haurwitz (lon, lat,u_wind, v_wind)
      CALL turnwi(u_wind,v_wind,u_tmp,v_tmp,lon_in,lat_in,lon,lat,0.0_r8,-0.5_r8*pi+rotation_angle,-1)
      fout = v_tmp
    case(11)
      fout = 1.0E-8_r8
    case(12)
      !
      ! Terminator chemistry initial condition
      !
      k1 = 1.0_r8*max(0.d0,sin(lat)*sin(k1_lat_center) + cos(lat)*cos(k1_lat_center)*cos(lon-k1_lon_center))
      k2 = 1._r8

      r = k1 / (4._r8*k2)
      det = sqrt(r*r + 2._r8*cly_constant*r)

      fout  = (det-r)
!      fout = cly_constant/2._r8 - (det-r)/2._r8
    case(13)
      !
      ! Terminator chemistry initial condition
      !
      k1 = 1.0_r8*max(0.d0,sin(lat)*sin(k1_lat_center) + cos(lat)*cos(k1_lat_center)*cos(lon-k1_lon_center))
      k2 = 1._r8

      r = k1 / (4._r8*k2)
      det = sqrt(r*r + 2._r8*cly_constant*r)

!      fout  = (det-r)
      fout = cly_constant/2._r8 - (det-r)/2._r8
    case default
!      call endrun("Illegal funcnum_arg in test_func")
!      fount = 1.0_r8
    end select
  end function test_func

  function test_wind(lat, lon, iwind) result(fout)
   use cam_abortutils, only: endrun
    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    integer,  intent(in) :: iwind

    real(r8)             :: fout


    fout = 0
  end function test_wind


  SUBROUTINE regrot(pxreg,pyreg,pxrot,pyrot,pxcen,pycen,kcall)
    use physconst, only: pi
!
!----------------------------------------------------------------------
!
!*    conversion between regular and rotated spherical coordinates.
!*
!*    pxreg     longitudes of the regular coordinates
!*    pyreg     latitudes of the regular coordinates
!*    pxrot     longitudes of the rotated coordinates
!*    pyrot     latitudes of the rotated coordinates
!*              all coordinates given in degrees n (negative for s)
!*              and degrees e (negative values for w)
!*    pxcen     regular longitude of the south pole of the rotated grid
!*    pycen     regular latitude of the south pole of the rotated grid
!*
!*    kcall=-1: find regular as functions of rotated coordinates.
!*    kcall= 1: find rotated as functions of regular coordinates.
!
!-----------------------------------------------------------------------
!
      integer kxdim,kydim,kx,ky,kcall
      real(r8) :: pxreg,pyreg,&
                  pxrot,pyrot,&
                  pxcen,pycen
!
!-----------------------------------------------------------------------
!
      real(r8) zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg,zcyreg, &
               zsyrot,zcyrot,zcxrot,zsxrot,zpi,zpih
      integer jy,jx

      zpih = pi*0.5_r8
!
      !----------------------------------------------------------------------
!
      zsycen = SIN((pycen+zpih))
      zcycen = COS((pycen+zpih))
!
      IF (kcall.eq.1) then
!
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc
         zsyrot = max(zsyrot,-1.0_r8)
         zsyrot = min(zsyrot,+1.0_r8)
         !
         pyrot = ASIN(zsyrot)
         !
         zcyrot = COS(pyrot)
         zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot
         zcxrot = max(zcxrot,-1.0_r8)
         zcxrot = min(zcxrot,+1.0_r8)
         zsxrot = zcyreg*zsxmxc/zcyrot
         !
         pxrot = ACOS(zcxrot)
         !
         IF (zsxrot < 0.0_r8) then
           pxrot = -pxrot
         end IF
               !
      ELSEIF (kcall.eq.-1) then
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         zsyreg = zcycen*zsyrot + zsycen*zcyrot*zcxrot
         zsyreg = max(zsyreg,-1.0_r8)
         zsyreg = min(zsyreg,+1.0_r8)
         !
         pyreg = ASIN(zsyreg)
         !
         zcyreg = COS(pyreg)
         zcxmxc = (zcycen*zcyrot*zcxrot -&
              zsycen*zsyrot)/zcyreg
         zcxmxc = max(zcxmxc,-1.0_r8)
         zcxmxc = min(zcxmxc,+1.0_r8)
         zsxmxc = zcyrot*zsxrot/zcyreg
         zxmxc  = ACOS(zcxmxc)
         IF (zsxmxc < 0.0_r8) then
           zxmxc = -zxmxc
         end IF
         !
         pxreg = zxmxc + pxcen
         !
      ELSE
         WRITE(6,'(1x,''invalid kcall in regrot'')')
         STOP
      ENDIF
    END SUBROUTINE regrot

  SUBROUTINE turnwi(puarg,pvarg,pures,pvres,pxreg,pyreg,pxrot,pyrot,pxcen,pycen,kcall)
    use physconst, only: pi
    !
    !-----------------------------------------------------------------------
    !
    !*    turn horizontal velocity components between regular and
    !*    rotated spherical coordinates.
    !
    !*    puarg : input u components
    !*    pvarg : input v components
    !*    pures : output u components
    !*    pvres : output v components
    !*    pa    : transformation coefficients
    !*    pb    :    -"-
    !*    pc    :    -"-
    !*    pd    :    -"-
    !*    pxreg : regular longitudes
    !*    pyreg : regular latitudes
    !*    pxrot : rotated longitudes
    !*    pyrot : rotated latitudes
    !*    kxdim              : dimension in the x (longitude) direction
    !*    kydim              : dimension in the y (latitude) direction
    !*    kx                 : number of gridpoints in the x direction
    !*    ky                 : number of gridpoints in the y direction
    !*    pxcen              : regular longitude of the south pole of the
    !*                         transformed grid
    !*    pycen              : regular latitude of the south pole of the
    !*                         transformed grid
    !*
    !*    kcall < 0          : find wind components in regular coordinates
    !*                         from wind components in rotated coordinates
    !*    kcall > 0          : find wind components in rotated coordinates
    !*                         from wind components in regular coordinates
    !*    note that all coordinates are given in degrees n and degrees e.
    !*       (negative values for s and w)
    !
    !-----------------------------------------------------------------------

    integer kxdim,kydim,kx,ky,kcall
    real(r8) puarg,pvarg,    &
         pures,pvres,    &
         pa,   pb,       &
         pc,   pd,       &
         pxreg,pyreg,    &
         pxrot,pyrot
      real(r8) pxcen,pycen
      !
      !-----------------------------------------------------------------------
      !
      integer jy,jx
      real(r8) zpih,zsyc,zcyc,zsxreg,zcxreg,zsyreg,zcyreg,zxmxc,&
           zsxmxc,zcxmxc,zsxrot,zcxrot,zsyrot,zcyrot
      !
      !-----------------------------------------------------------------------
      !
      IF (kcall.eq.1) then
        zpih = pi*0.5_r8
        zsyc = SIN(pycen+zpih)
        zcyc = COS(pycen+zpih)
        !
        zsxreg = SIN(pxreg)
        zcxreg = COS(pxreg)
        zsyreg = SIN(pyreg)
        zcyreg = COS(pyreg)
        !
        zxmxc  = pxreg - pxcen
        zsxmxc = SIN(zxmxc)
        zcxmxc = COS(zxmxc)
        !
        zsxrot = SIN(pxrot)
        zcxrot = COS(pxrot)
        zsyrot = SIN(pyrot)
        zcyrot = COS(pyrot)
        !
        pa = zcyc*zsxmxc*zsxrot + zcxmxc*zcxrot
        pb = zcyc*zcxmxc*zsyreg*zsxrot - zsyc*zcyreg*zsxrot - &
             zsxmxc*zsyreg*zcxrot
        pc = zsyc*zsxmxc/zcyrot
        pd = (zsyc*zcxmxc*zsyreg + zcyc*zcyreg)/zcyrot
        !
        pures = pa*puarg + pb*pvarg
        pvres = pc*puarg + pd*pvarg
      ELSEIF (kcall.eq.-1) then
        zpih = pi*0.5_r8
        zsyc = SIN(pycen+zpih)
        zcyc = COS(pycen+zpih)
        !
        zsxreg = SIN(pxreg)
        zcxreg = COS(pxreg)
        zsyreg = SIN(pyreg)
        zcyreg = COS(pyreg)
        !
        zxmxc  = pxreg - pxcen
        zsxmxc = SIN(zxmxc)
        zcxmxc = COS(zxmxc)
        !
        zsxrot = SIN(pxrot)
        zcxrot = COS(pxrot)
        zsyrot = SIN(pyrot)
        zcyrot = COS(pyrot)
        !
        pa = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot
        pb = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot -&
             zcxmxc*zsxrot*zsyrot
        pc =-zsyc*zsxrot/zcyreg
        pd = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg
        !
        pures = pa*puarg + pb*pvarg
        pvres = pc*puarg + pd*pvarg
      ELSE
        write(6,'(1x,''invalid kcall in turnwi'')')
        STOP
      ENDIF
    END SUBROUTINE turnwi

  SUBROUTINE Rossby_Haurwitz (lon, lat,u_wind, v_wind)
    use physconst, only: rearth
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(in)  :: lon,                                     & ! longitude in radians
                               lat                                        ! latitude in radians
                                                                          ! both coefficients 'a' and 'b' are needed at the full model level
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind                                     ! meridional wind in m/s

!-----------------------------------------------------------------------
!     test case parameters
!-----------------------------------------------------------------------
      real(r8),parameter :: u0      = 50._r8,                          &   ! reference wind
                            n       = 4._r8                                ! wavenumber

!-----------------------------------------------------------------------
!     local
!-----------------------------------------------------------------------
      real(r8) :: tmp1, tmp2, tmp3, KK, MM
      real(r8) :: sin_lat, cos_lat, sin_slat, cos_slat

!-----------------------------------------------------------------------
!     initialize the wind components
!-----------------------------------------------------------------------
      MM      = u0/(n*rearth)   ! parameter M
      KK      = u0/(n*rearth)   ! parameter K


      cos_lat = cos(lat)
      sin_lat = sin(lat)
      tmp1 = rearth * MM * cos_lat
      tmp2 = rearth * KK * cos_lat**(n-1._r8)*(n*sin_lat**2 - cos_lat**2)
      tmp3 = -rearth * KK * n * cos_lat**(n-1._r8) * sin_lat
      u_wind = tmp1 + tmp2 * cos(n*lon)
      v_wind = tmp3 * sin(n*lon)
  end subroutine Rossby_Haurwitz

#endif
end module test_fvm_mapping
