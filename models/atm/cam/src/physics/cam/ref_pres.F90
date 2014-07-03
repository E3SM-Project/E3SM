module ref_pres
!--------------------------------------------------------------------------
! 
! Provides access to reference pressures for use by the physics
! parameterizations.  The pressures are provided by the dynamical core
! since it currently determines the grid used by the physics.
! 
! Note that the init method for this module is called before the init
! method in physpkg; therefore, most physics modules can use these
! reference pressures during their init phases.
! 
!--------------------------------------------------------------------------

use shr_kind_mod, only: r8=>shr_kind_r8
use ppgrid,       only: pver, pverp

implicit none
public
save

! Reference pressures (Pa)
real(r8), protected :: pref_edge(pverp)     ! Layer edges
real(r8), protected :: pref_mid(pver)       ! Layer midpoints
real(r8), protected :: pref_mid_norm(pver)  ! Layer midpoints normalized by
                                            ! surface pressure ('eta' coordinate)

real(r8), protected :: ptop_ref             ! Top of model
real(r8), protected :: psurf_ref            ! Surface pressure

! Number of top levels using pure pressure representation
integer, protected :: num_pr_lev

! Pressure used to set troposphere cloud physics top (Pa)
real(r8), protected :: trop_cloud_top_press = 0._r8
! Top level for troposphere cloud physics
integer, protected :: trop_cloud_top_lev

! Pressure used to set MAM process top (Pa)
real(r8), protected :: clim_modal_aero_top_press = 0._r8
! Top level for MAM processes that impact climate
integer, protected :: clim_modal_aero_top_lev

! Molecular diffusion is calculated only if the model top is below this
! pressure (Pa).
real(r8), protected :: do_molec_press = 0.1_r8
! Pressure used to set bottom of molecular diffusion region (Pa).
real(r8), protected :: molec_diff_bot_press = 50._r8
! Flag for molecular diffusion, and molecular diffusion level indices.
logical, protected :: do_molec_diff = .false.
integer, protected :: ntop_molec = 1
integer, protected :: nbot_molec = 0

!====================================================================================
contains
!====================================================================================

subroutine ref_pres_readnl(nlfile)

   use spmd_utils,      only: masterproc
   use abortutils,      only: endrun
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'ref_pres_readnl'

   namelist /ref_pres_nl/ trop_cloud_top_press, clim_modal_aero_top_press,&
        do_molec_press, molec_diff_bot_press
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'ref_pres_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, ref_pres_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! Check that top for modal aerosols is not lower than
      ! top for clouds.
      if (clim_modal_aero_top_press > trop_cloud_top_press) &
           call endrun("ERROR: clim_modal_aero_top press must be less &
           &than or equal to trop_cloud_top_press.")
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(trop_cloud_top_press,            1 , mpir8,   0, mpicom)
   call mpibcast(clim_modal_aero_top_press,       1 , mpir8,   0, mpicom)
   call mpibcast(do_molec_press,                  1 , mpir8,   0, mpicom)
   call mpibcast(molec_diff_bot_press,            1 , mpir8,   0, mpicom)
#endif

end subroutine ref_pres_readnl

!====================================================================================

subroutine ref_pres_init

  use dyn_grid,     only: dyn_grid_get_pref

  ! Get reference pressures from the dynamical core.
  call dyn_grid_get_pref(pref_edge, pref_mid, num_pr_lev)

  ptop_ref = pref_edge(1)
  psurf_ref = pref_edge(pverp)

  pref_mid_norm = pref_mid/psurf_ref

  ! Find level corresponding to the top of troposphere clouds.
  trop_cloud_top_lev = press_lim_idx(trop_cloud_top_press, &
       top=.true.)

  ! Find level corresponding to the top for MAM processes.
  clim_modal_aero_top_lev = press_lim_idx(clim_modal_aero_top_press, &
       top=.true.)

  ! Find level corresponding to the molecular diffusion bottom.
  do_molec_diff = (ptop_ref < do_molec_press)
  if (do_molec_diff) then
     nbot_molec = press_lim_idx(molec_diff_bot_press, &
          top=.false.)
  end if

end subroutine ref_pres_init

!====================================================================================

! Convert pressure limiters to the appropriate level.
pure function press_lim_idx(p, top) result(k_lim)
  ! Pressure
  real(r8), intent(in) :: p
  ! Is this a top or bottom limit?
  logical,  intent(in) :: top
  integer :: k_lim, k

  if (top) then
     k_lim = pver+1
     do k = 1, pver
        if (pref_mid(k) > p) then
           k_lim = k
           exit
        end if
     end do
  else
     k_lim = 0
     do k = pver, 1, -1
        if (pref_mid(k) < p) then
           k_lim = k
           exit
        end if
     end do
  end if

end function press_lim_idx

!====================================================================================

end module ref_pres
