module ref_pres
!------------------------------------------------------------------------------------
! 
! Provides access to reference pressures for use by the physics
! parameterizations.  The pressures are provided by the dynamical core
! since it currently determines the grid used by the physics.
! 
!------------------------------------------------------------------------------------

use shr_kind_mod, only: r8=>shr_kind_r8
use ppgrid,       only: pver, pverp
use dyn_grid,     only: dyn_grid_get_pref

implicit none
public
save

real(r8) :: pref_edge(pverp)     ! reference pressure at layer edges (Pa)
real(r8) :: pref_mid(pver)       ! reference pressure at layer midpoints (Pa)
real(r8) :: pref_mid_norm(pver)  ! reference pressure at layer midpoints normalized
                                 ! by the surface pressure ('eta' coordinate)

real(r8) :: ptop_ref             ! reference pressure at top of model (Pa)
real(r8) :: psurf_ref            ! reference surface pressure (Pa)

integer  :: num_pr_lev           ! number of top levels using pure pressure representation

! Pressure used to set troposphere cloud physics top (Pa)
real(r8) :: trop_cloud_top_press = 0._r8
! Top level for troposphere cloud physics
integer  :: trop_cloud_top_lev

! Pressure used to set MAM process top (Pa)
real(r8) :: clim_modal_aero_top_press = 0._r8
! Top level for MAM processes that impact climate
integer  :: clim_modal_aero_top_lev

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

   namelist /ref_pres_nl/ trop_cloud_top_press, clim_modal_aero_top_press
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
#endif

end subroutine ref_pres_readnl

!====================================================================================

subroutine ref_pres_init

   call dyn_grid_get_pref(pref_edge, pref_mid, num_pr_lev)

   ptop_ref = pref_edge(1)
   psurf_ref = pref_edge(pverp)

   pref_mid_norm = pref_mid/psurf_ref

   ! Find level corresponding to the top of troposphere clouds.
   trop_cloud_top_lev = minloc(abs(pref_mid - trop_cloud_top_press), 1)

   ! Find level corresponding to the top for MAM processes.
   clim_modal_aero_top_lev = minloc(abs(pref_mid - clim_modal_aero_top_press), 1)

end subroutine ref_pres_init

!====================================================================================

end module ref_pres
