module dadadj_cam

! CAM interfaces for the dry adiabatic adjustment parameterization

use shr_kind_mod,    only: r8=>shr_kind_r8, cs=>shr_kind_cs
use ppgrid,          only: pcols, pver, pverp
use constituents,    only: pcnst
use physconst,       only: cappav, cpairv, pi
use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
use phys_control,    only: q3d_is_on
use cam_abortutils,  only: endrun
use cam_logfile,     only: iulog
use error_messages,  only: handle_errmsg

use spmd_utils,      only: masterproc, masterprocid, mpicom, mpi_integer
use namelist_utils,  only: find_group_name
use units,           only: getunit, freeunit

use dadadj,          only: dadadj_initial, dadadj_calc

implicit none
private
save

public :: &
   dadadj_readnl, &
   dadadj_init, &
   dadadj_tend

! Namelist variables
integer :: dadadj_nlvdry = 3  ! number of layers from top of model to apply the adjustment
integer :: dadadj_niter = 15  ! number of iterations for convergence

!===============================================================================
contains
!===============================================================================

subroutine dadadj_readnl(filein)

   !!!character(len=cs), intent(in) :: filein ! Input namelist filename
   character(len=*), intent(in) :: filein ! Input namelist filename

   namelist /dadadj_nl/ dadadj_nlvdry, dadadj_niter

   integer :: unitn, ierr
   character(len=*), parameter :: sub='dadadj_readnl'
   !------------------------------------------------------------------

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
        write(*,*) '### unitn = ',unitn,' filein = ',trim(filein)
      open(unitn, file=trim(filein), status='old')
      call find_group_name(unitn, 'dadadj_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, dadadj_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun( sub//':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(dadadj_nlvdry, 1, mpi_integer, masterprocid, mpicom)
   call mpibcast(dadadj_niter, 1, mpi_integer, masterprocid, mpicom)
#endif

   call dadadj_initial(dadadj_nlvdry, dadadj_niter)

   !if (masterproc .and. .not. use_simple_phys) then
   if (masterproc .and. .not. q3d_is_on) then
      write(iulog,*)'Dry adiabatic adjustment applied to top N layers; N=', &
                    dadadj_nlvdry
      write(iulog,*)'Dry adiabatic adjustment number of iterations for convergence =', &
                    dadadj_niter
   end if

end subroutine dadadj_readnl


!===============================================================================

subroutine dadadj_init()
    use cam_history,   only: addfld

    call addfld('DADADJ_PD', (/ 'lev' /), 'A', 'probability', 'dry adiabatic adjustment probability')

end subroutine dadadj_init


!===============================================================================

subroutine dadadj_tend(dt, state, ptend)
   use cam_history,   only: outfld

   real(r8),                  intent(in)  :: dt         ! Time step [s]
   type(physics_state),       intent(in)  :: state      ! Physics state variables
   type(physics_ptend),       intent(out) :: ptend      ! parameterization tendencies

   logical :: lq(pcnst)
   real(r8) :: dadpdf(pcols, pver)
   integer :: ncol, lchnk, icol_err
   character(len=128) :: errstring  ! Error string

    ncol  = state%ncol
    lchnk = state%lchnk
    lq(:) = .FALSE.
    lq(1) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'dadadj', ls=.true., lq=lq)

    ! use the ptend components for temporary storate and copy state info for input to
    ! dadadj_calc which directly updates the temperature and moisture input arrays.

    ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
    ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

    call dadadj_calc( &
       ncol, state%pmid, state%pint, state%pdel, cappav(:,:,lchnk), ptend%s, &
       ptend%q(:,:,1), dadpdf, icol_err)
       
    call outfld('DADADJ_PD',  dadpdf(:ncol,:),  ncol, lchnk)

    if (icol_err > 0) then
       ! error exit
       write(errstring, *) &
          'dadadj_calc: No convergence in column at lat,lon:', &
          state%lat(icol_err)*180._r8/pi, state%lon(icol_err)*180._r8/pi
       call handle_errmsg(errstring, subname="dadadj_tend")
    end if

    ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/dt * cpairv(:ncol,:,lchnk)
    ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/dt

end subroutine dadadj_tend

!===============================================================================
end module dadadj_cam
