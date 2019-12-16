module forcing_interface

  use iso_c_binding,  only: c_int, c_bool, c_double, c_ptr, c_f_pointer
  use dimensions_mod, only: nlev, nlevp, np
  use element_mod,    only: element_t
  use kinds,          only: real_kind 
  use hybvcoord_mod,  only: hvcoord_t 
  use parallel_mod,   only: abortmp

  implicit none

  type(hvcoord_t) :: hvcoord

  type(element_t), allocatable :: elem(:)
  real (kind=real_kind), pointer, dimension(:,:,:,:)     :: ps, ft, fvtheta, fphi
  real (kind=real_kind), pointer, dimension(:,:,:,:,:)   :: q, fq, w, vtheta, dp, phinh, fm
  real (kind=real_kind), pointer, dimension(:,:,:,:,:,:) :: qdp, v

  public :: init_f90
  public :: set_forcing_pointers_f90
  public :: tracers_forcing_f90
  public :: cleanup_f90

contains

  subroutine init_f90 (num_elems, hyai, hybi, hyam, hybm, gradphis, ps0, qsz) bind(c)
    use dimensions_mod, only: nelemd, qsize
    use element_state,  only: allocate_element_arrays, setup_element_pointers_ie
    !
    ! Inputs
    !
    real (kind=real_kind), intent(in) :: hyai(nlevp), hybi(nlevp)
    real (kind=real_kind), intent(in) :: hyam(nlev), hybm(nlev)
    real (kind=real_kind), intent(in) :: gradphis(np,np,2,num_elems)
    real (kind=real_kind), intent(in) :: ps0
    integer (kind=c_int),  intent(in) :: num_elems, qsz
    !
    ! Locals
    !
    integer :: ie

    hvcoord%hyai = hyai
    hvcoord%hybi = hybi
    hvcoord%hyam = hyam
    hvcoord%hybm = hybm
    hvcoord%ps0 = ps0
    qsize = qsz

    nelemd = num_elems

    call allocate_element_arrays(num_elems)

    allocate (elem(num_elems))
    do ie=1,num_elems
      call setup_element_pointers_ie(ie, elem(ie)%state, elem(ie)%derived, elem(ie)%accum)
      elem(ie)%derived%gradphis = gradphis(:,:,:,ie)
    enddo

  end subroutine init_f90

  subroutine set_forcing_pointers_f90 (q_ptr, fq_ptr, qdp_ptr,    &
                                       v_ptr, w_ptr, vtheta_ptr,  &
                                       dp_ptr, phinh_ptr, ps_ptr, &
                                       fm_ptr, ft_ptr,            &
                                       fvtheta_ptr, fphi_ptr) bind(c)
    use dimensions_mod, only: nelemd, qsize_d
    use element_state,  only: timelevels
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: q_ptr, fq_ptr, qdp_ptr
    type (c_ptr), intent(in) :: v_ptr, w_ptr, vtheta_ptr, dp_ptr, phinh_ptr, ps_ptr
    type (c_ptr), intent(in) :: fm_ptr, ft_ptr, fvtheta_ptr, fphi_ptr
    !
    ! Locals
    !
    integer :: ie

    call c_f_pointer (q_ptr,      q,       [np,np,nlev,qsize_d,nelemd])
    call c_f_pointer (fq_ptr,     fq,      [np,np,nlev,qsize_d,nelemd])
    call c_f_pointer (qdp_ptr,    qdp,     [np,np,nlev,qsize_d,2,nelemd])
    call c_f_pointer (v_ptr,      v,       [np,np,2,nlev, timelevels,nelemd])
    call c_f_pointer (w_ptr,      w,       [np,np,  nlevp,timelevels,nelemd])
    call c_f_pointer (vtheta_ptr, vtheta,  [np,np,  nlev, timelevels,nelemd])
    call c_f_pointer (dp_ptr,     dp,      [np,np,  nlev, timelevels,nelemd])
    call c_f_pointer (phinh_ptr,  phinh,   [np,np,  nlevp,timelevels,nelemd])
    call c_f_pointer (ps_ptr,     ps,      [np,np,        timelevels,nelemd])
    call c_f_pointer (fm_ptr,     fm,      [np,np,3,nlev, nelemd])
    call c_f_pointer (ft_ptr,     ft,      [np,np,  nlev, nelemd])
    call c_f_pointer (fvtheta_ptr,fvtheta, [np,np,  nlev, nelemd])
    call c_f_pointer (fphi_ptr,   fphi,    [np,np,  nlevp,nelemd])
  end subroutine set_forcing_pointers_f90

  subroutine dynamics_forcing_f90(dt,np1) bind(c)
    use dimensions_mod,   only: nelemd
    use prim_advance_mod, only: applyCAMforcing_dynamics
    !
    ! Inputs
    !
    real (kind=c_double),  intent(in) :: dt
    integer (kind=c_int),  intent(in) :: np1
    !
    ! Locals
    !
    integer :: ie

    do ie=1,nelemd
      ! Copy inputs from C ptrs
      elem(ie)%derived%FM      = fm(:,:,:,:,ie)
      elem(ie)%derived%FT      = ft(:,:,:,ie)
      elem(ie)%derived%FVTheta = fvtheta(:,:,:,ie)
      elem(ie)%derived%FPHI    = fphi(:,:,:,ie)
    enddo

    call applyCAMforcing_dynamics(elem,hvcoord,np1,dt,1,nelemd)

    do ie=1,nelemd
      ! Copy outputs to C ptrs
      v(:,:,:,:,:,ie)    = elem(ie)%state%v        
      w(:,:,:,:,ie)      = elem(ie)%state%w_i      
      vtheta(:,:,:,:,ie) = elem(ie)%state%vtheta_dp
      phinh(:,:,:,:,ie)  = elem(ie)%state%phinh_i  
    enddo
  end subroutine dynamics_forcing_f90

  subroutine tracers_forcing_f90(dt,np1,np1_qdp,hydrostatic,moist) bind(c)
    use control_mod,      only: use_moisture, theta_hydrostatic_mode
    use dimensions_mod,   only: nelemd
    use prim_driver_base, only: applyCAMforcing_tracers
    !
    ! Inputs
    !
    real (kind=c_double),  intent(in) :: dt
    integer (kind=c_int),  intent(in) :: np1,np1_qdp
    logical (kind=c_bool), intent(in) :: hydrostatic, moist
    !
    ! Locals
    !
    integer :: ie

    theta_hydrostatic_mode = hydrostatic
    use_moisture = moist

    do ie=1,nelemd
      ! Copy inputs from C ptrs
      elem(ie)%state%q    = q(:,:,:,:,ie)
      elem(ie)%derived%FQ = fq(:,:,:,:,ie)
      elem(ie)%state%qdp  = qdp(:,:,:,:,:,ie)

      elem(ie)%state%v         = v(:,:,:,:,:,ie)
      elem(ie)%state%w_i       = w(:,:,:,:,ie)
      elem(ie)%state%dp3d      = dp(:,:,:,:,ie)
      elem(ie)%state%vtheta_dp = vtheta(:,:,:,:,ie)
      elem(ie)%state%phinh_i   = phinh(:,:,:,:,ie)
      elem(ie)%state%ps_v      = ps(:,:,:,ie)

      elem(ie)%derived%FM      = fm(:,:,:,:,ie)
      elem(ie)%derived%FT      = ft(:,:,:,ie)
      elem(ie)%derived%FVTheta = fvtheta(:,:,:,ie)
      elem(ie)%derived%FPHI    = fphi(:,:,:,ie)

      ! Apply forcing
      call applyCAMforcing_tracers(elem(ie),hvcoord,np1,np1_qdp,dt,.false.)

      ! Copy outputs to C ptrs
      q(:,:,:,:,ie)     = elem(ie)%state%q    
      fq(:,:,:,:,ie)    = elem(ie)%derived%FQ 
      qdp(:,:,:,:,:,ie) = elem(ie)%state%qdp  

      v(:,:,:,:,:,ie)    = elem(ie)%state%v        
      w(:,:,:,:,ie)      = elem(ie)%state%w_i      
      dp(:,:,:,:,ie)     = elem(ie)%state%dp3d     
      vtheta(:,:,:,:,ie) = elem(ie)%state%vtheta_dp
      phinh(:,:,:,:,ie)  = elem(ie)%state%phinh_i  
      ps(:,:,:,ie)       = elem(ie)%state%ps_v     
                           
      fm(:,:,:,:,ie)     = elem(ie)%derived%FM     
      ft(:,:,:,ie)       = elem(ie)%derived%FT     
      fvtheta(:,:,:,ie)  = elem(ie)%derived%FVTheta
      fphi(:,:,:,ie)     = elem(ie)%derived%FPHI   
    enddo
  end subroutine tracers_forcing_f90

  subroutine cleanup_f90 () bind(c)
    use element_state, only: deallocate_element_arrays

    call deallocate_element_arrays()
    deallocate(elem)
  end subroutine cleanup_f90

end module forcing_interface
