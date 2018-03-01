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
  use shr_kind_mod, only : cxx => SHR_KIND_CXX
  use glc_elevclass_mod, only : glc_get_num_elevation_classes, glc_get_elevation_class, &
       glc_elevclass_as_string, glc_all_elevclass_strings, GLC_ELEVCLASS_STRLEN, &
       GLC_ELEVCLASS_ERR_NONE, GLC_ELEVCLASS_ERR_TOO_LOW, &
       GLC_ELEVCLASS_ERR_TOO_HIGH, glc_errcode_to_string
  use mct_mod
  use seq_map_type_mod, only : seq_map
  use seq_map_mod, only : seq_map_map
  use shr_sys_mod, only : shr_sys_abort

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: map_lnd2glc  ! map one field from LND -> GLC grid

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: get_glc_elevation_classes ! get the elevation class of each point on the glc grid
  private :: map_bare_land             ! remap the field of interest for the bare land "elevation class"
  private :: map_ice_covered           ! remap the field of interest for all elevation classes (excluding bare land)

contains

  !-----------------------------------------------------------------------
  subroutine map_lnd2glc(l2x_l, landfrac_l, g2x_g, fieldname, &
       mapper, l2x_g)
    !
    ! !DESCRIPTION:
    ! Maps one field from the LND grid to the GLC grid.
    !
    ! Mapping is done with a multiplication by landfrac on the source grid, with
    ! normalization.
    !
    ! Sets the given field within l2x_g, leaving the rest of l2x_g untouched.
    !
    ! Assumes that l2x_l contains fields like:
    ! - fieldname00
    ! - fieldname01
    ! - fieldname02
    ! - etc.
    !
    ! and also:
    ! - Sl_topo00
    ! - Sl_topo01
    ! - Sl_topo02
    ! - etc.
    !
    ! and l2x_g contains a field named 'fieldname'
    !
    ! Assumes that landfrac_l contains the field:
    ! - lfrac: land fraction on the land grid
    !
    ! Assumes that g2x_g contains the following fields:
    ! - Sg_ice_covered: whether each glc grid cell is ice-covered (0 or 1)
    ! - Sg_topo: ice topographic height on the glc grid
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect)  , intent(in)    :: l2x_l      ! lnd -> cpl fields on the land grid
    type(mct_aVect)  , intent(in)    :: landfrac_l ! lfrac field on the land grid
    type(mct_aVect)  , intent(in)    :: g2x_g      ! glc -> cpl fields on the glc grid
    character(len=*) , intent(in)    :: fieldname  ! name of the field to map
    type(seq_map)    , intent(inout) :: mapper
    type(mct_aVect)  , intent(inout) :: l2x_g      ! lnd -> cpl fields on the glc grid
    !
    ! !LOCAL VARIABLES:

    ! fieldname with trailing blanks removed
    character(len=:), allocatable :: fieldname_trimmed

    ! number of points on the GLC grid
    integer :: lsize_g

    ! data for bare land on the GLC grid
    ! needs to be a pointer to satisfy the MCT interface
    real(r8), pointer :: data_g_bareland(:)

    ! data for ice-covered regions on the GLC grid
    ! needs to be a pointer to satisfy the MCT interface
    real(r8), pointer :: data_g_ice_covered(:)

    ! final data on the GLC grid
    ! needs to be a pointer to satisfy the MCT interface
    real(r8), pointer :: data_g(:)

    ! whether each point on the glc grid is ice-covered (1) or ice-free (0)
    ! needs to be a pointer to satisfy the MCT interface
    real(r8), pointer :: glc_ice_covered(:)

    ! ice topographic height on the glc grid
    ! needs to be a pointer to satisfy the MCT interface
    real(r8), pointer :: glc_topo(:)

    ! elevation class on the glc grid
    ! 0 implies bare ground (no ice)
    integer, allocatable :: glc_elevclass(:)

    character(len=*), parameter :: subname = 'map_lnd2glc'
    !-----------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Initialize temporary arrays and other local variables
    ! ------------------------------------------------------------------------

    lsize_g = mct_aVect_lsize(l2x_g)

    allocate(data_g_ice_covered(lsize_g))
    allocate(data_g_bareland(lsize_g))
    allocate(data_g(lsize_g))

    fieldname_trimmed = trim(fieldname)

    ! ------------------------------------------------------------------------
    ! Extract necessary fields from g2x_g
    ! ------------------------------------------------------------------------

    allocate(glc_ice_covered(lsize_g))
    allocate(glc_topo(lsize_g))
    call mct_aVect_exportRattr(g2x_g, 'Sg_ice_covered', glc_ice_covered)
    call mct_aVect_exportRattr(g2x_g, 'Sg_topo', glc_topo)

    ! ------------------------------------------------------------------------
    ! Determine elevation class of each glc point
    ! ------------------------------------------------------------------------

    allocate(glc_elevclass(lsize_g))
    call get_glc_elevation_classes(glc_ice_covered, glc_topo, glc_elevclass)

    ! ------------------------------------------------------------------------
    ! Map elevation class 0 (bare land) and ice elevation classes
    ! ------------------------------------------------------------------------

    call map_bare_land(l2x_l, landfrac_l, fieldname_trimmed, mapper, data_g_bareland)

    ! Start by setting the output data equal to the bare land value everywhere; this will
    ! later get overwritten in places where we have ice
    !
    ! TODO(wjs, 2015-01-20) This implies that we pass data to CISM even in places that
    ! CISM says is ocean (so CISM will ignore the incoming value). This differs from the
    ! current glint implementation, which sets acab and artm to 0 over ocean (although
    ! notes that this could lead to a loss of conservation). Figure out how to handle
    ! this case.
    data_g(:) = data_g_bareland(:)

    ! Map the SMB to ice-covered cells
    call map_ice_covered(l2x_l, landfrac_l, fieldname_trimmed, &
         glc_topo, mapper, data_g_ice_covered)

    where (glc_elevclass /= 0)
       data_g = data_g_ice_covered
    end where

    ! ------------------------------------------------------------------------
    ! Set field in output attribute vector
    ! ------------------------------------------------------------------------

    call mct_aVect_importRattr(l2x_g, fieldname_trimmed, data_g)

    ! ------------------------------------------------------------------------
    ! Clean up
    ! ------------------------------------------------------------------------

    deallocate(data_g_ice_covered)
    deallocate(data_g_bareland)
    deallocate(data_g)
    deallocate(glc_ice_covered)
    deallocate(glc_topo)
    deallocate(glc_elevclass)

  end subroutine map_lnd2glc

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

  !-----------------------------------------------------------------------
  subroutine map_bare_land(l2x_l, landfrac_l, fieldname, mapper, data_g_bare_land)
    !
    ! !DESCRIPTION:
    ! Remaps the field of interest for the bare land "elevation class".
    !
    ! Puts the output in data_g_bare_land, which should already be allocated to have size
    ! equal to the number of GLC points that this processor is responsible for.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect)  , intent(in)    :: l2x_l      ! lnd -> cpl fields on the land grid
    type(mct_aVect)  , intent(in)    :: landfrac_l ! lfrac field on the land grid
    character(len=*) , intent(in)    :: fieldname  ! name of the field to map (should have NO trailing blanks)
    type(seq_map)    , intent(inout) :: mapper
    real(r8), pointer, intent(inout) :: data_g_bare_land(:)
    !
    ! !LOCAL VARIABLES:
    character(len=:), allocatable :: elevclass_as_string
    character(len=:), allocatable :: fieldname_bare_land
    integer :: lsize_g  ! number of points for attribute vectors on the glc grid
    type(mct_aVect) :: l2x_g_bare_land  ! temporary attribute vector holding the remapped field for bare land

    character(len=*), parameter :: subname = 'map_bare_land'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(associated(data_g_bare_land), __FILE__, __LINE__)

    lsize_g = size(data_g_bare_land)
    elevclass_as_string = glc_elevclass_as_string(0)
    fieldname_bare_land = fieldname // elevclass_as_string
    call mct_aVect_init(l2x_g_bare_land, rList = fieldname_bare_land, lsize = lsize_g)

    call seq_map_map(mapper = mapper, av_s = l2x_l, av_d = l2x_g_bare_land, &
         fldlist = fieldname_bare_land, &
         norm = .true., &
         avwts_s = landfrac_l, &
         avwtsfld_s = 'lfrac')
    call mct_aVect_exportRattr(l2x_g_bare_land, fieldname_bare_land, data_g_bare_land)

    call mct_aVect_clean(l2x_g_bare_land)

  end subroutine map_bare_land

  !-----------------------------------------------------------------------
  subroutine map_ice_covered(l2x_l, landfrac_l, fieldname, &
       topo_g, mapper, data_g_ice_covered)

    !
    ! !DESCRIPTION:
    ! Remaps the field of interest from the land grid (in multiple elevation classes)
    ! to the glc grid
    !
    ! Puts the output in data_g_ice_covered, which should already be allocated to have size
    ! equal to the number of GLC points that this processor is responsible for.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect)  , intent(in)    :: l2x_l      ! lnd -> cpl fields on the land grid
    type(mct_aVect)  , intent(in)    :: landfrac_l ! lfrac field on the land grid
    character(len=*) , intent(in)    :: fieldname  ! name of the field to map (should have NO trailing blanks)
    real(r8)         , intent(in)    :: topo_g(:)  ! topographic height for each point on the glc grid
    type(seq_map)    , intent(inout) :: mapper
    real(r8)         , intent(out)   :: data_g_ice_covered(:)  ! field remapped to glc grid

    ! !LOCAL VARIABLES:

    character(len=*), parameter :: toponame = 'Sl_topo'  ! base name for topo fields in l2x_l;
    ! actual names will have elevation class suffix

    character(len=GLC_ELEVCLASS_STRLEN), allocatable :: all_elevclass_strings(:)
    character(len=:), allocatable :: elevclass_as_string
    character(len=:), allocatable :: fieldname_ec
    character(len=:), allocatable :: toponame_ec
    character(len=:), allocatable :: fieldnamelist
    character(len=:), allocatable :: toponamelist
    character(len=:), allocatable :: totalfieldlist

    integer :: nEC           ! number of elevation classes
    integer :: lsize_g       ! number of cells on glc grid
    integer :: n, ec
    integer :: strlen

    real(r8) :: elev_l, elev_u  ! lower and upper elevations in interpolation range
    real(r8) :: d_elev          ! elev_u - elev_l

    type(mct_aVect) :: l2x_g_temp  ! temporary attribute vector holding the remapped fields for this elevation class

    real(r8), pointer :: tmp_field_g(:)  ! must be a pointer to satisfy the MCT interface
    real, pointer :: data_g_EC(:,:)    ! remapped field in each glc cell, in each EC
    real, pointer :: topo_g_EC(:,:)    ! remapped topo in each glc cell, in each EC

    ! 1 is probably enough, but use 10 to be safe, in case the length of the delimiter
    ! changes
    integer, parameter :: extra_len_for_list_merge = 10

    character(len=*), parameter :: subname = 'map_ice_covered'
    !-----------------------------------------------------------------------

    lsize_g = size(data_g_ice_covered)
    nEC = glc_get_num_elevation_classes()
    SHR_ASSERT_FL((size(topo_g) == lsize_g), __FILE__, __LINE__)

    ! ------------------------------------------------------------------------
    ! Create temporary vectors
    ! ------------------------------------------------------------------------

    allocate(tmp_field_g(lsize_g))
    allocate(data_g_EC  (lsize_g,nEC))
    allocate(topo_g_EC  (lsize_g,nEC))

    ! ------------------------------------------------------------------------
    ! Make a string that concatenates all EC levels of field, as well as the topo
    ! The resulting list will look something like this:
    !    'Flgl_qice01:Flgl_qice02: ... :Flgl_qice10:Sl_topo01:Sl_topo02: ... :Sltopo10'
    ! ------------------------------------------------------------------------

    allocate(all_elevclass_strings(1:glc_get_num_elevation_classes()))
    all_elevclass_strings = glc_all_elevclass_strings(include_zero = .false.)
    fieldnamelist = shr_string_listFromSuffixes( &
         suffixes = all_elevclass_strings, &
         strBase  = fieldname)
    toponamelist = shr_string_listFromSuffixes( &
         suffixes = all_elevclass_strings, &
         strBase  = toponame)
    strlen = len_trim(fieldnamelist) + len_trim(toponamelist) + extra_len_for_list_merge
    allocate(character(len=strlen) :: totalfieldlist)
    call shr_string_listMerge(fieldnamelist, toponamelist, totalfieldlist )

    ! ------------------------------------------------------------------------
    ! Make a temporary attribute vector.
    ! For each grid cell on the land grid, this attribute vector contains the field and
    ! topo values for all ECs.
    ! ------------------------------------------------------------------------
    call mct_aVect_init(l2x_g_temp, rList = totalfieldlist, lsize = lsize_g)

    ! ------------------------------------------------------------------------
    ! Remap all these fields from the land (source) grid to the glc (destination) grid.
    ! ------------------------------------------------------------------------

    call seq_map_map(mapper = mapper, &
         av_s = l2x_l, &
         av_d = l2x_g_temp, &
         fldlist = totalfieldlist, &
         norm = .true., &
         avwts_s = landfrac_l, &
         avwtsfld_s = 'lfrac')

    ! ------------------------------------------------------------------------
    ! Export all elevation classes out of attribute vector and into local 2D arrays (xy,z)
    ! ------------------------------------------------------------------------

    do ec = 1, nEC
       elevclass_as_string = glc_elevclass_as_string(ec)
       fieldname_ec = fieldname // elevclass_as_string
       toponame_ec = toponame // elevclass_as_string
       call mct_aVect_exportRattr(l2x_g_temp, fieldname_ec, tmp_field_g)
       data_g_EC(:,ec) = real(tmp_field_g)
       call mct_aVect_exportRattr(l2x_g_temp, toponame_ec, tmp_field_g)
       topo_g_EC(:,ec) = real(tmp_field_g)
    enddo

    ! ------------------------------------------------------------------------
    ! Perform vertical interpolation of data onto ice sheet topography
    ! ------------------------------------------------------------------------

    data_g_ice_covered(:) = 0._r8

    do n = 1, lsize_g

       ! For each ice sheet point, find bounding EC values...
       if (topo_g(n) < topo_g_EC(n,1)) then
          ! lower than lowest mean EC elevation value
          data_g_ice_covered(n) = data_g_EC(n,1)

       else if (topo_g(n) >= topo_g_EC(n,nEC)) then
          ! higher than highest mean EC elevation value
          data_g_ice_covered(n) = data_g_EC(n,nEC)

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
                   data_g_ice_covered(n) = data_g_EC(n,ec-1) * 0.5_r8 &
                        + data_g_EC(n,ec)   * 0.5_r8
                else
                   data_g_ice_covered(n) = data_g_EC(n,ec-1) * (elev_u - topo_g(n)) / d_elev  &
                        + data_g_EC(n,ec)   * (topo_g(n) - elev_l) / d_elev
                end if

                exit
             end if
          end do
       end if  ! topo_g(n)
    end do  ! lsize_g

    ! ------------------------------------------------------------------------
    ! Clean up
    ! ------------------------------------------------------------------------

    deallocate(tmp_field_g)
    deallocate(data_g_EC)
    deallocate(topo_g_EC)

    call mct_aVect_clean(l2x_g_temp)

  end subroutine map_ice_covered

end module map_lnd2glc_mod
