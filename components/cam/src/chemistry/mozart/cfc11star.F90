!---------------------------------------------------------------------------------
! Manages the CFC11* for radiation 
!  4 Dec 2009 -- Francis Vitt created
!---------------------------------------------------------------------------------
module cfc11star

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  
  use physics_buffer, only : pbuf_add_field, dtype_r8
  use cam_abortutils,   only : endrun
  use ppgrid,       only : pcols, pver, begchunk, endchunk
  use spmd_utils,   only : masterproc

  implicit none
  save 

  private
  public :: register_cfc11star
  public :: update_cfc11star
  public :: init_cfc11star

  logical :: do_cfc11star
  character(len=16), parameter :: pbufname = 'CFC11STAR'
  integer :: pbf_idx
 
  integer, pointer :: cfc11_ndx
  integer, pointer :: cfc113_ndx
  integer, pointer :: ccl4_ndx
  integer, pointer :: ch3ccl3_ndx
  integer, pointer :: hcfc22_ndx
  integer, pointer :: cf2clbr_ndx
  integer, pointer :: cf3br_ndx
  integer, target :: indices(7)
  
  real(r8) :: rel_rf(7)

contains

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
  subroutine register_cfc11star

    implicit none
    
    call pbuf_add_field(pbufname,'global',dtype_r8,(/pcols,pver/),pbf_idx)

  endsubroutine register_cfc11star

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
  subroutine init_cfc11star(pbuf2d)
    use constituents, only : cnst_get_ind
    use cam_history,  only : addfld
    use infnan,       only : nan, assignment(=)
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field

    implicit none

    real(r8) :: real_nan

    real(r8), parameter :: cfc_rf(7)  =  (/ 0.25_r8, 0.30_r8, 0.13_r8, 0.06_r8, 0.20_r8, 0.30_r8, 0.32_r8 /)    ! W/m2/ppb
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    real_nan = nan

    cfc11_ndx   => indices(1)
    cfc113_ndx  => indices(2)
    ccl4_ndx    => indices(3)
    ch3ccl3_ndx => indices(4)
    hcfc22_ndx  => indices(5)
    cf2clbr_ndx => indices(6)
    cf3br_ndx   => indices(7)
    
    call cnst_get_ind('CFC11',  cfc11_ndx,   abort=.false.)
    call cnst_get_ind('CFC113', cfc113_ndx,  abort=.false.)
    call cnst_get_ind('CCL4',   ccl4_ndx,    abort=.false.)
    call cnst_get_ind('CH3CCL3',ch3ccl3_ndx, abort=.false.)
    call cnst_get_ind('HCFC22', hcfc22_ndx,  abort=.false.)
    call cnst_get_ind('CF2CLBR',cf2clbr_ndx, abort=.false.)
    call cnst_get_ind('CF3BR',  cf3br_ndx,   abort=.false.)

    do_cfc11star = all(indices(:)>0)

    if (.not.do_cfc11star) return

    call pbuf_set_field(pbuf2d, pbf_idx, real_nan)

    rel_rf(:) = cfc_rf(:) / cfc_rf(1)
    call addfld(pbufname,(/ 'lev' /),'A','kg/kg','cfc11star for radiation' )
    
    if (masterproc) then
       write(iulog,*) 'init_cfc11star: CFC11STAR is added to pbuf2d for radiation'
    endif
  end subroutine init_cfc11star

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
  subroutine update_cfc11star( pbuf2d, phys_state )
    use cam_history,  only : outfld
    use physics_types,only : physics_state
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

    implicit none

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


    integer :: lchnk, ncol
    integer :: c
    real(r8), pointer :: cf11star(:,:)

    if (.not.do_cfc11star) return
    
    do c = begchunk,endchunk
       lchnk = phys_state(c)%lchnk
       ncol = phys_state(c)%ncol

       call pbuf_get_field(pbuf_get_chunk(pbuf2d, lchnk), pbf_idx, cf11star)

       cf11star(:ncol,:) = &
            phys_state(c)%q(:ncol,:,cfc11_ndx)  * rel_rf(1) + &
            phys_state(c)%q(:ncol,:,cfc113_ndx) * rel_rf(2) + &
            phys_state(c)%q(:ncol,:,ccl4_ndx)   * rel_rf(3) + &
            phys_state(c)%q(:ncol,:,ch3ccl3_ndx)* rel_rf(4) + &
            phys_state(c)%q(:ncol,:,hcfc22_ndx) * rel_rf(5) + &
            phys_state(c)%q(:ncol,:,cf2clbr_ndx)* rel_rf(6) + &
            phys_state(c)%q(:ncol,:,cf3br_ndx)  * rel_rf(7)

       call outfld( pbufname, cf11star(:ncol,:), ncol, lchnk )

    enddo

  endsubroutine update_cfc11star

end module cfc11star
