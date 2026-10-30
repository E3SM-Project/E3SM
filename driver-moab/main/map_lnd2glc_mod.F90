module map_lnd2glc_mod

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module contains routines for mapping fields from the LND grid (separated by GLC
  ! elevation class) onto the GLC grid
  !
  ! For high-level design, see:
  ! https://docs.google.com/document/d/1H_SuK6SfCv1x6dK91q80dFInPbLYcOkUj_iAa6WRnqQ/edit

#include "shr_assert.h"
  use seq_comm_mct, only: CPLID, GLCID, logunit
  use shr_kind_mod, only : r8 => shr_kind_r8
  use glc_elevclass_mod, only : glc_get_elevation_class, &
       glc_elevclass_as_string, glc_all_elevclass_strings, GLC_ELEVCLASS_STRLEN, &
       GLC_ELEVCLASS_ERR_NONE, GLC_ELEVCLASS_ERR_TOO_LOW, &
       GLC_ELEVCLASS_ERR_TOO_HIGH, glc_errcode_to_string
  use mct_mod
  use shr_sys_mod, only : shr_sys_abort

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  ! array-based pieces of the algorithm, used by the moab path in prep_glc_mod
  ! (the horizontal map is done there with one batched seq_map_map call on tags;
  ! these routines provide the elevation-class assignment and the vertical
  ! interpolation on plain arrays fetched from the glc mesh tags)
  public :: get_glc_elevation_classes  ! get the elevation class of each point on the glc grid
  public :: map_lnd2glc_vertical_interp ! vertically interpolate per-EC data to the ice sheet topography

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------


contains

  !-----------------------------------------------------------------------
  subroutine map_lnd2glc_vertical_interp(topo_g, topo_g_EC, data_g_EC, data_g_bareland, &
       glc_elevclass, data_g)
    !
    ! !DESCRIPTION:
    ! Vertically interpolate a field, already horizontally mapped to the glc grid in
    ! each elevation class, onto the actual ice sheet topography. This is the
    ! array-based equivalent of the combination of map_bare_land / map_ice_covered
    ! (minus the horizontal maps, which the caller has already done): the output is
    ! the bare-land (EC 0) value where the glc cell is ice-free, and the linear
    ! vertical interpolation between bounding elevation classes where ice-covered.
    !
    ! All arrays are on the glc grid decomposition. data_g_EC and topo_g_EC hold the
    ! horizontally-mapped field and Sl_topo for elevation classes 1..nEC.
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: topo_g(:)          ! ice topographic height on the glc grid
    real(r8), intent(in)  :: topo_g_EC(:,:)     ! mapped per-EC topo (lsize_g, nEC)
    real(r8), intent(in)  :: data_g_EC(:,:)     ! mapped per-EC field (lsize_g, nEC)
    real(r8), intent(in)  :: data_g_bareland(:) ! mapped bare-land (EC 0) field
    integer , intent(in)  :: glc_elevclass(:)   ! elevation class of each glc point (0 = bare)
    real(r8), intent(out) :: data_g(:)          ! result on the glc grid
    !
    ! !LOCAL VARIABLES:
    integer :: lsize_g       ! number of cells on glc grid
    integer :: nEC           ! number of elevation classes
    integer :: n, ec
    real(r8) :: elev_l, elev_u  ! lower and upper elevations in interpolation range
    real(r8) :: d_elev          ! elev_u - elev_l

    character(len=*), parameter :: subname = 'map_lnd2glc_vertical_interp'
    !-----------------------------------------------------------------------

    lsize_g = size(data_g)
    nEC = size(data_g_EC, 2)
    SHR_ASSERT_FL((size(topo_g) == lsize_g), __FILE__, __LINE__)
    SHR_ASSERT_FL((size(glc_elevclass) == lsize_g), __FILE__, __LINE__)

    do n = 1, lsize_g

       if (glc_elevclass(n) == 0) then
          ! bare land (ice-free) point: use the bare-land value
          data_g(n) = data_g_bareland(n)

       ! For each ice sheet point, find bounding EC values...
       else if (topo_g(n) < topo_g_EC(n,1)) then
          ! lower than lowest mean EC elevation value
          data_g(n) = data_g_EC(n,1)

       else if (topo_g(n) >= topo_g_EC(n,nEC)) then
          ! higher than highest mean EC elevation value
          data_g(n) = data_g_EC(n,nEC)

       else
          ! do linear interpolation of data in the vertical
          do ec = 2, nEC
             if (topo_g(n) < topo_g_EC(n, ec)) then
                elev_l = topo_g_EC(n, ec-1)
                elev_u = topo_g_EC(n, ec)
                d_elev = elev_u - elev_l
                if (d_elev <= 0) then
                   ! This shouldn't happen, but handle it in case it does. In this case,
                   ! let's arbitrarily use the mean of the two elevation classes, rather
                   ! than the weighted mean.
                   write(logunit,*) subname//' WARNING: topo diff between elevation classes <= 0'
                   write(logunit,*) 'n, ec, elev_l, elev_u = ', n, ec, elev_l, elev_u
                   write(logunit,*) 'Simply using mean of the two elevation classes,'
                   write(logunit,*) 'rather than the weighted mean.'
                   data_g(n) = data_g_EC(n,ec-1) * 0.5_r8 &
                        + data_g_EC(n,ec)   * 0.5_r8
                else
                   data_g(n) = data_g_EC(n,ec-1) * (elev_u - topo_g(n)) / d_elev  &
                        + data_g_EC(n,ec)   * (topo_g(n) - elev_l) / d_elev
                end if

                exit
             end if
          end do
       end if  ! topo_g(n)
    end do  ! lsize_g

  end subroutine map_lnd2glc_vertical_interp

  !-----------------------------------------------------------------------
  subroutine get_glc_elevation_classes(glc_ice_covered, glc_topo, glc_elevclass)
    !
    ! !DESCRIPTION:
    ! Get the elevation class of each point on the glc grid.
    !
    ! For grid cells that are ice-free, the elevation class is set to 0.
    !
    ! All arguments (glc_ice_covered, glc_topo and glc_elevclass) must be the same size.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: glc_ice_covered(:) ! ice-covered (1) vs. ice-free (0)
    real(r8), intent(in)  :: glc_topo(:)        ! ice topographic height
    integer , intent(out) :: glc_elevclass(:)   ! elevation class
    !
    ! !LOCAL VARIABLES:
    integer :: npts
    integer :: glc_pt
    integer :: err_code

    ! Tolerance for checking whether ice_covered is 0 or 1
    real(r8), parameter :: ice_covered_tol = 1.e-13

    character(len=*), parameter :: subname = 'get_glc_elevation_classes'
    !-----------------------------------------------------------------------

    npts = size(glc_elevclass)
    SHR_ASSERT_FL((size(glc_ice_covered) == npts), __FILE__, __LINE__)
    SHR_ASSERT_FL((size(glc_topo) == npts), __FILE__, __LINE__)

    do glc_pt = 1, npts
       if (abs(glc_ice_covered(glc_pt) - 1._r8) < ice_covered_tol) then
          ! This is an ice-covered point

          call glc_get_elevation_class(glc_topo(glc_pt), glc_elevclass(glc_pt), err_code)
          if ( err_code == GLC_ELEVCLASS_ERR_NONE .or. &
               err_code == GLC_ELEVCLASS_ERR_TOO_LOW .or. &
               err_code == GLC_ELEVCLASS_ERR_TOO_HIGH) then
             ! These are all acceptable "errors" - it is even okay for these purposes if
             ! the elevation is lower than the lower bound of elevation class 1, or
             ! higher than the upper bound of the top elevation class.

             ! Do nothing
          else
             write(logunit,*) subname, ': ERROR getting elevation class for ', glc_pt
             write(logunit,*) glc_errcode_to_string(err_code)
             call shr_sys_abort(subname//': ERROR getting elevation class')
          end if
       else if (abs(glc_ice_covered(glc_pt) - 0._r8) < ice_covered_tol) then
          ! This is a bare land point (no ice)
          glc_elevclass(glc_pt) = 0
       else
          ! glc_ice_covered is some value other than 0 or 1
          ! The lnd -> glc downscaling code would need to be reworked if we wanted to
          ! handle a continuous fraction between 0 and 1.
          write(logunit,*) subname, ': ERROR: glc_ice_covered must be 0 or 1'
          write(logunit,*) 'glc_pt, glc_ice_covered = ', glc_pt, glc_ice_covered(glc_pt)
          call shr_sys_abort(subname//': ERROR: glc_ice_covered must be 0 or 1')
       end if
    end do

  end subroutine get_glc_elevation_classes

end module map_lnd2glc_mod
