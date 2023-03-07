module microp_driver

!-------------------------------------------------------------------------------------------------------
!
! Driver for CAM microphysics parameterizations
!
!-------------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use ppgrid,         only: pver
use physics_types,  only: physics_state, physics_ptend, physics_tend,  &
                          physics_ptend_copy, physics_ptend_sum
use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
use phys_control,   only: phys_getopts

use micro_mg_cam,   only: micro_mg_cam_readnl, micro_mg_cam_register, &
                          micro_mg_cam_implements_cnst, micro_mg_cam_init_cnst, &
                          micro_mg_cam_init, micro_mg_cam_tend
use micro_p3_interface, only: micro_p3_init, micro_p3_register, micro_p3_tend, &
                              micro_p3_implements_cnst, micro_p3_init_cnst, &
                              micro_p3_readnl
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use perf_mod,       only: t_startf, t_stopf

implicit none
private
save

public :: &
   microp_driver_readnl,          &
   microp_driver_register,        &
   microp_driver_init_cnst,       &
   microp_driver_implements_cnst, &
   microp_driver_init,            &
   microp_driver_tend

character(len=16)  :: microp_scheme   ! Microphysics scheme
!===============================================================================
contains
!===============================================================================

subroutine microp_driver_readnl(nlfile)

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Read in namelist for microphysics scheme
   !-----------------------------------------------------------------------

   call phys_getopts(microp_scheme_out=microp_scheme)

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_readnl(nlfile)
   case ('P3')
      call micro_p3_readnl(nlfile)
   ! microp_driver doesn't handle these other options
   case ('RK')
      continue
   case ('off')
      continue
   case default
      call endrun('microp_driver_readnl:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_readnl

subroutine microp_driver_register

   ! Register microphysics constituents and fields in the physics buffer.
   !-----------------------------------------------------------------------


   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_register()
   case ('P3')
      call micro_p3_register()
   ! microp_driver doesn't handle these other options
   case ('RK')
      continue
   case ('off')
      continue
   case default
      call endrun('microp_driver_register:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_register

!===============================================================================

function microp_driver_implements_cnst(name)

   ! Return true if specified constituent is implemented by the
   ! microphysics package

   character(len=*), intent(in) :: name        ! constituent name
   logical :: microp_driver_implements_cnst    ! return value

   ! Local workspace
   integer :: m
   !-----------------------------------------------------------------------

   microp_driver_implements_cnst = .false.

   select case (microp_scheme)
   case ('MG')
      microp_driver_implements_cnst = micro_mg_cam_implements_cnst(name)
   case ('P3')
      microp_driver_implements_cnst = micro_p3_implements_cnst(name)
   ! microp_driver doesn't handle these other options
   case ('RK')
      continue
   case ('off')
      continue
   case default
      call endrun('microp_driver_implements_cnst:: unrecognized microp_scheme')
   end select

end function microp_driver_implements_cnst

!===============================================================================

subroutine microp_driver_init_cnst(name, q, gcid)

   ! Initialize the microphysics constituents, if they are
   ! not read from the initial file.

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_init_cnst(name, q, gcid)
   case ('P3')
      call micro_p3_init_cnst(name, q)
   ! microp_driver doesn't handle these other options
   case ('RK')
      continue
   case ('off')
      continue
   case default
      call endrun('microp_driver_init_cnst:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_init_cnst

!===============================================================================

subroutine microp_driver_init(pbuf2d)

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! Initialize the microphysics parameterizations
   !-----------------------------------------------------------------------

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_init(pbuf2d)
   case ('P3')
      call micro_p3_init(pbuf2d)
   ! microp_driver doesn't handle these other options
   case ('RK')
      continue
   case ('off')
      continue
   case default
      call endrun('microp_driver_init:: unrecognized microp_scheme')
   end select


end subroutine microp_driver_init

!===============================================================================

subroutine microp_driver_tend(state, ptend, dtime, pbuf)

   ! Call the microphysics parameterization run methods.

   ! Input arguments

   type(physics_state), intent(in)    :: state       ! State variables
   type(physics_ptend), intent(out)   :: ptend       ! Package tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in)  :: dtime                    ! Timestep

   ! Local variables

   integer :: lchnk
   integer :: ncol

   !======================================================================

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Call MG Microphysics

   select case (microp_scheme)
   case ('MG')
      call t_startf('microp_mg_cam_tend')
      call micro_mg_cam_tend(state, ptend, dtime, pbuf)
      call t_stopf('microp_mg_cam_tend')
   case ('P3')
      call t_startf('microp_p3_tend')
      call micro_p3_tend(state, ptend, dtime, pbuf)
      call t_stopf('microp_p3_tend')
   ! microp_driver doesn't handle these other options
   case ('RK')
      continue
   case ('off')
      continue
   case default
      call endrun('microp_driver_tend:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_tend

end module microp_driver
