module map_lnd2glc_mod

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module contains routines for mapping fields from the LND grid (separated by GLC
  ! elevation class) onto the GLC grid
  !
  ! For high-level design, see:
  ! https://docs.google.com/document/d/1sjsaiPYsPJ9A7dVGJIHGg4rVIY2qF5aRXbNzSXVAafU/edit?usp=sharing

#include "shr_assert.h"
  use seq_comm_mct, only : logunit
  use shr_kind_mod, only : r8 => shr_kind_r8
  use glc_elevclass_mod, only : glc_get_num_elevation_classes, glc_get_elevation_class, &
       glc_elevclass_as_string, GLC_ELEVCLASS_ERR_NONE, GLC_ELEVCLASS_ERR_TOO_LOW, &
       GLC_ELEVCLASS_ERR_TOO_HIGH, glc_errcode_to_string
  use mct_mod
  use seq_map_type_mod, only : seq_map
  use seq_map_mod, only : seq_map_map
  use vertical_gradient_calculator_base, only : vertical_gradient_calculator_base_type
  use shr_log_mod, only : errMsg => shr_log_errMsg
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
  private :: map_one_elevation_class   ! remap the field of interest for one ice elevation class

contains

  !-----------------------------------------------------------------------
  subroutine map_lnd2glc(l2x_l, landfrac_l, g2x_g, fieldname, gradient_calculator, &
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
    class(vertical_gradient_calculator_base_type), intent(inout) :: gradient_calculator
    type(seq_map)    , intent(inout) :: mapper
    type(mct_aVect)  , intent(inout) :: l2x_g      ! lnd -> cpl fields on the glc grid
    !
    ! !LOCAL VARIABLES:

    ! fieldname with trailing blanks removed
    character(len=:), allocatable :: fieldname_trimmed

    ! index for looping over elevation classes
    integer :: elevclass

    ! number of points on the GLC grid
    integer :: lsize_g

    ! data for one elevation class on the GLC grid
    ! needs to be a pointer to satisfy the MCT interface
    real(r8), pointer :: data_g_oneEC(:)

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
    allocate(data_g_oneEC(lsize_g))
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
    ! Map elevation class 0 (bare land)
    ! ------------------------------------------------------------------------

    call map_bare_land(l2x_l, landfrac_l, fieldname_trimmed, mapper, data_g_oneEC)

    ! Start by setting the output data equal to the bare land value everywhere; this will
    ! later get overwritten in places where we have ice
    !
    ! TODO(wjs, 2015-01-20) This implies that we pass data to CISM even in places that
    ! CISM says is ocean (so CISM will ignore the incoming value). This differs from the
    ! current glint implementation, which sets acab and artm to 0 over ocean (although
    ! notes that this could lead to a loss of conservation). Figure out how to handle
    ! this case.
    data_g(:) = data_g_oneEC(:)

    ! ------------------------------------------------------------------------
    ! Map ice elevation classes
    ! ------------------------------------------------------------------------

    call gradient_calculator%calc_gradients()
    do elevclass = 1, glc_get_num_elevation_classes()
       call map_one_elevation_class(l2x_l, landfrac_l, fieldname_trimmed, elevclass, &
            gradient_calculator, glc_topo, mapper, data_g_oneEC)

       where (glc_elevclass == elevclass)
          data_g = data_g_oneEC
       end where
    end do

    ! ------------------------------------------------------------------------
    ! Set field in output attribute vector
    ! ------------------------------------------------------------------------

    call mct_aVect_importRattr(l2x_g, fieldname_trimmed, data_g)

    ! ------------------------------------------------------------------------
    ! Clean up
    ! ------------------------------------------------------------------------

    deallocate(data_g_oneEC)
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
    SHR_ASSERT((size(glc_ice_covered) == npts), errMsg(__FILE__, __LINE__))
    SHR_ASSERT((size(glc_topo) == npts), errMsg(__FILE__, __LINE__))

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

    SHR_ASSERT(associated(data_g_bare_land), errMsg(__FILE__, __LINE__))

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
  subroutine map_one_elevation_class(l2x_l, landfrac_l, fieldname, elevclass, &
       gradient_calculator, topo_g, mapper, data_g_thisEC)
    !
    ! !DESCRIPTION:
    ! Remaps the field of interest for a single ice elevation class.
    !
    ! Puts the output in data_g_thisEC, which should already be allocated to have size
    ! equal to the number of GLC points that this processor is responsible for.
    !
    ! To do this remapping, we remap the field adjusted by the vertical gradient. That is,
    ! rather than mapping data_l itself, we instead remap:
    !
    !   data_l + (vertical_gradient_l) * (topo_g - topo_l)
    !
    ! (where _l denotes quantities on the land grid, _g on the glc grid)
    !
    ! However, in order to do the remapping with existing routines, we do this by
    ! performing two separate remappings:
    !
    !   (1) Remap (data_l - vertical_gradient_l * topo_l); put result in partial_remap_g
    !       (note: in variables in the code, the parenthesized term is called partial_remap,
    !       either on the land or glc grid)
    !
    !   (2) Remap vertical_gradient_l; put result in vertical_gradient_g
    !
    ! Then data_g = partial_remap_g + topo_g * vertical_gradient_g
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect)  , intent(in)    :: l2x_l      ! lnd -> cpl fields on the land grid
    type(mct_aVect)  , intent(in)    :: landfrac_l ! lfrac field on the land grid
    character(len=*) , intent(in)    :: fieldname  ! name of the field to map (should have NO trailing blanks)
    integer          , intent(in)    :: elevclass  ! elevation class index to map
    class(vertical_gradient_calculator_base_type), intent(in) :: gradient_calculator
    real(r8)         , intent(in)    :: topo_g(:)  ! topographic height for each point on the glc grid
    type(seq_map)    , intent(inout) :: mapper
    real(r8)         , intent(out)   :: data_g_thisEC(:)
    !
    ! !LOCAL VARIABLES:

    ! Fields contained in the temporary attribute vectors:
    character(len=*), parameter :: partial_remap_tag = 'partial_remap'
    character(len=*), parameter :: vertical_gradient_tag = 'vertical_gradient'
    character(len=*), parameter :: attr_tags = &
         partial_remap_tag // ':' // vertical_gradient_tag

    ! Base name for the topo fields in l2x_l. Actual fields will have an elevation class suffix.
    character(len=*), parameter :: toponame = 'Sl_topo'

    character(len=:), allocatable :: elevclass_as_string
    character(len=:), allocatable :: fieldname_ec
    character(len=:), allocatable :: toponame_ec
    integer :: lsize_l  ! number of points for attribute vectors on the land grid
    integer :: lsize_g  ! number of points for attribute vectors on the glc grid
    type(mct_aVect) :: l2x_l_temp
    type(mct_aVect) :: l2x_g_temp  ! temporary attribute vector holding the remapped fields for this elevation class

    ! Note that arrays passed to MCT routines need to be pointers
    ! Temporary fields on the land grid:
    real(r8), pointer :: data_l(:)
    real(r8), pointer :: topo_l(:)
    real(r8), pointer :: vertical_gradient_l(:)
    real(r8), pointer :: partial_remap_l(:)
    ! Temporary fields on the glc grid:
    real(r8), pointer :: vertical_gradient_g(:)
    real(r8), pointer :: partial_remap_g(:)

    character(len=*), parameter :: subname = 'map_one_elevation_class'
    !-----------------------------------------------------------------------

    lsize_g = size(data_g_thisEC)
    SHR_ASSERT((size(topo_g) == lsize_g), errMsg(__FILE__, __LINE__))

    ! ------------------------------------------------------------------------
    ! Create temporary attribute vectors
    ! ------------------------------------------------------------------------

    lsize_l = mct_aVect_lsize(l2x_l)
    call mct_aVect_init(l2x_l_temp, rList = attr_tags, lsize = lsize_l)
    call mct_aVect_init(l2x_g_temp, rList = attr_tags, lsize = lsize_g)

    ! ------------------------------------------------------------------------
    ! Create fields to remap on the source (land) grid
    ! ------------------------------------------------------------------------

    allocate(data_l(lsize_l))
    allocate(topo_l(lsize_l))
    allocate(vertical_gradient_l(lsize_l))
    allocate(partial_remap_l(lsize_l))

    elevclass_as_string = glc_elevclass_as_string(elevclass)
    fieldname_ec = fieldname // elevclass_as_string
    toponame_ec = toponame // elevclass_as_string

    call mct_aVect_exportRattr(l2x_l, fieldname_ec, data_l)
    call mct_aVect_exportRattr(l2x_l, toponame_ec, topo_l)
    call gradient_calculator%get_gradients_one_class(elevclass, vertical_gradient_l)
    partial_remap_l = data_l - (vertical_gradient_l * topo_l)

    call mct_aVect_importRattr(l2x_l_temp, partial_remap_tag, partial_remap_l)
    call mct_aVect_importRattr(l2x_l_temp, vertical_gradient_tag, vertical_gradient_l)

    ! ------------------------------------------------------------------------
    ! Remap to destination (glc) grid
    ! ------------------------------------------------------------------------

    call seq_map_map(mapper = mapper, av_s = l2x_l_temp, av_d = l2x_g_temp, &
         norm = .true., &
         avwts_s = landfrac_l, &
         avwtsfld_s = 'lfrac')

    ! ------------------------------------------------------------------------
    ! Compute final field on the destination (glc) grid
    ! ------------------------------------------------------------------------

    allocate(partial_remap_g(lsize_g))
    allocate(vertical_gradient_g(lsize_g))

    call mct_aVect_exportRattr(l2x_g_temp, partial_remap_tag, partial_remap_g)
    call mct_aVect_exportRattr(l2x_g_temp, vertical_gradient_tag, vertical_gradient_g)

    data_g_thisEC = partial_remap_g + (topo_g * vertical_gradient_g)

    ! ------------------------------------------------------------------------
    ! Clean up
    ! ------------------------------------------------------------------------

    call mct_aVect_clean(l2x_g_temp)
    call mct_aVect_clean(l2x_l_temp)
    deallocate(data_l)
    deallocate(topo_l)
    deallocate(vertical_gradient_l)
    deallocate(partial_remap_l)
    deallocate(vertical_gradient_g)
    deallocate(partial_remap_g)

  end subroutine map_one_elevation_class



end module map_lnd2glc_mod
