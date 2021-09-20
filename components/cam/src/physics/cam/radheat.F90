
module radheat
!-----------------------------------------------------------------------
!
! Purpose:  Provide an interface to convert shortwave and longwave
!           radiative heating terms into net heating.
!
!           This module provides a hook to allow incorporating additional
!           radiative terms (eUV heating and nonLTE longwave cooling).
! 
! Original version: B.A. Boville
!-----------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use ppgrid,        only: pcols, pver
use physics_types, only: physics_state, physics_ptend, physics_ptend_init

use physics_buffer, only : physics_buffer_desc

implicit none
private
save

! Public interfaces
public  &
   radheat_readnl,        &!
   radheat_init,          &!
   radheat_timestep_init, &!
   radheat_tend            ! return net radiative heating

!===============================================================================
contains
!===============================================================================

subroutine radheat_readnl(nlfile)

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! No options for this version of radheat; this is just a stub.

end subroutine radheat_readnl

!================================================================================================

subroutine radheat_init(pref_mid)

   use pmgrid, only: plev
   use physics_buffer, only : physics_buffer_desc

   real(r8), intent(in) :: pref_mid(plev)


end subroutine radheat_init

!================================================================================================

subroutine radheat_timestep_init (state, pbuf2d)
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


end subroutine radheat_timestep_init

!================================================================================================

subroutine radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
                        fsnt, flns, flnt, asdir, net_flx)
#if ( defined OFFLINE_DYN )
   use metdata, only: met_rlx, met_srf_feedback
#endif
!++BEH
    use prescribed_radheat, only: has_presc_radheat
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_index
!--BEH
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!-----------------------------------------------------------------------

! Arguments
   type(physics_state), intent(in)  :: state             ! Physics state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_ptend), intent(out) :: ptend             ! indivdual parameterization tendencie
!++BEH
!B   real(r8),            intent(in)  :: qrl(pcols,pver)   ! longwave heating
!B   real(r8),            intent(in)  :: qrs(pcols,pver)   ! shortwave heating
!B   real(r8),            intent(in)  :: fsns(pcols)       ! Surface solar absorbed flux
!B   real(r8),            intent(in)  :: fsnt(pcols)       ! Net column abs solar flux at model top
!B   real(r8),            intent(in)  :: flns(pcols)       ! Srf longwave cooling (up-down) flux
!B   real(r8),            intent(in)  :: flnt(pcols)       ! Net outgoing lw flux at model top
   real(r8),            intent(inout) :: qrl(pcols,pver)   ! longwave heating
   real(r8),            intent(inout) :: qrs(pcols,pver)   ! shortwave heating
   real(r8),            intent(inout) :: fsns(pcols)       ! Surface solar absorbed flux
   real(r8),            intent(inout) :: fsnt(pcols)       ! Net column abs solar flux at model top
   real(r8),            intent(inout) :: flns(pcols)       ! Srf longwave cooling (up-down) flux
   real(r8),            intent(inout) :: flnt(pcols)       ! Net outgoing lw flux at model top
!--BEH
   real(r8),            intent(in)  :: asdir(pcols)      ! shortwave, direct albedo
   real(r8),            intent(out) :: net_flx(pcols)  


! Local variables
   integer :: i, k
   integer :: ncol
!++BEH
   integer :: index_qrs, index_qrl
   integer :: index_fsnt, index_flnt, index_fsns, index_flns
   real(r8), pointer, dimension(:,:) :: qrs2d, qrl2d
   real(r8), pointer, dimension(:,:) :: fsnt2d, flnt2d, fsns2d, flns2d
!--BEH
!-----------------------------------------------------------------------

   ncol = state%ncol

   call physics_ptend_init(ptend,state%psetcols, 'cam_radheat', ls=.true.)

!++BEH
   if ( has_presc_radheat ) then
      index_qrs   = pbuf_get_index('p_QRS')
      index_qrl   = pbuf_get_index('p_QRL')
      index_fsnt  = pbuf_get_index('p_FSNT')
      index_flnt  = pbuf_get_index('p_FLNT')
      index_fsns  = pbuf_get_index('p_FSNS')
      index_flns  = pbuf_get_index('p_FLNS')

      call pbuf_get_field(pbuf, index_qrs,  qrs2d)
      call pbuf_get_field(pbuf, index_qrl,  qrl2d)
      call pbuf_get_field(pbuf, index_fsnt, fsnt2d)
      call pbuf_get_field(pbuf, index_flnt, flnt2d)
      call pbuf_get_field(pbuf, index_fsns, fsns2d)
      call pbuf_get_field(pbuf, index_flns, flns2d)
      do i = 1, ncol
         fsnt(i) = fsnt2d(i,1)
         flnt(i) = flnt2d(i,1)
         fsns(i) = fsns2d(i,1)
         flns(i) = flns2d(i,1)
         do k = 1, pver
            qrs(i,k) = qrs2d(i,k)
            qrl(i,k) = qrl2d(i,k)
         end do
      end do
   end if
!--BEH

#if ( defined OFFLINE_DYN )
   ptend%s(:ncol,:) = 0._r8
   do k = 1,pver
     if (met_rlx(k) < 1._r8 .or. met_srf_feedback) then
       ptend%s(:ncol,k) = (qrs(:ncol,k) + qrl(:ncol,k))
     endif
   enddo 
#else
   ptend%s(:ncol,:) = (qrs(:ncol,:) + qrl(:ncol,:))
#endif

   do i = 1, ncol
      net_flx(i) = fsnt(i) - fsns(i) - flnt(i) + flns(i)
   end do

end subroutine radheat_tend
end module radheat
