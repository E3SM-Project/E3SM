module map_glc2lnd_mod

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module contains routines for mapping fields from the GLC grid onto the LND grid
  ! (separated by GLC elevation class)
  !
  ! For high-level design, see:
  ! https://docs.google.com/document/d/1sjsaiPYsPJ9A7dVGJIHGg4rVIY2qF5aRXbNzSXVAafU/edit?usp=sharing

#include "shr_assert.h"
  use seq_comm_mct, only : logunit
  use shr_kind_mod, only : r8 => shr_kind_r8
  use glc_elevclass_mod, only : glc_get_num_elevation_classes, glc_get_elevation_class, &
       glc_mean_elevation_virtual, glc_elevclass_as_string, &
       GLC_ELEVCLASS_ERR_NONE, GLC_ELEVCLASS_ERR_TOO_LOW, &
       GLC_ELEVCLASS_ERR_TOO_HIGH, glc_errcode_to_string
  use mct_mod
  use seq_map_type_mod, only : seq_map
  use seq_map_mod, only : seq_map_map
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use shr_sys_mod, only : shr_sys_abort

  
  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: map_glc2lnd_ec  ! map all fields from GLC -> LND grid that need to be separated by elevation class

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: get_glc_elevation_classes     ! get elevation class of each glc cell
  private :: get_frac_this_ec              ! get fraction in a given elevation class
  private :: set_topo_in_virtual_columns
  private :: make_aVect_frac_times_icemask
  
  character(len=*), parameter :: frac_times_icemask_field = 'Sg_frac_times_icemask'
  
contains

  !-----------------------------------------------------------------------
  subroutine map_glc2lnd_ec(g2x_g, &
       frac_field, topo_field, icemask_field, extra_fields, &
       mapper, g2x_l)
    !
    ! !DESCRIPTION:
    ! Maps fields from the GLC grid to the LND grid that need to be separated by
    ! elevation class.
    !
    ! Maps frac_field, topo_field, plus all fields defined in extra_fields. extra_fields
    ! should be a colon-delimited list of fields, giving the field name in the g2x_g
    ! attribute vector (i.e., without the elevation class suffixes).
    !
    ! Assumes that g2x_g contains:
    ! - frac_field
    ! - topo_field
    ! - icemask_field (Note: this is NOT mapped here, but is needed as an input to the mapping)
    ! - each field in extra_fields
    !
    ! Assumes that g2x_l contains:
    ! - <frac_field>00, <frac_field>01, <frac_field>02, ...
    ! - <topo_field>00, <topo_field>01, <topo_field>02, ...
    ! - And similarly for each field in extra_fields
    !
    ! Currently assumes that all fields are mapped using the same mapper, which should be
    ! a conservative mapper (i.e., a flux mapper).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect), intent(in) :: g2x_g
    character(len=*), intent(in) :: frac_field    ! name of field in g2x_g containing glc ice fraction
    character(len=*), intent(in) :: topo_field    ! name of field in g2x_g containing glc topo
    character(len=*), intent(in) :: icemask_field ! name of field in g2x_g containing ice mask
    character(len=*), intent(in) :: extra_fields
    type(seq_map), intent(inout) :: mapper
    type(mct_aVect), intent(inout) :: g2x_l
    
    !
    ! !LOCAL VARIABLES:
    integer :: lsize_g
    integer :: lsize_l

    ! The following need to be pointers to satisfy the MCT interface:
    real(r8), pointer :: glc_frac(:)  ! total ice fraction in each glc cell
    real(r8), pointer :: glc_topo(:)  ! topographic height of each glc cell
    real(r8), pointer :: glc_frac_this_ec(:)  ! ice fraction in this elevation class, for eachglc cell

    integer , allocatable :: glc_elevclass(:)  ! elevation class of each glc cell (assuming cell is ice-covered)
    integer :: n
    character(len=:), allocatable :: elevclass_as_string
    character(len=:), allocatable :: frac_field_ec  ! field name: frac_field with elev class suffix
    character(len=len(extra_fields)+100) :: fields_to_map
    character(len=2*len(extra_fields)+100) :: fields_to_map_ec  ! fields_to_map with elev class suffixes
    integer :: num_fields_to_map
    
    ! attribute vector holding glc fraction in one elev class, on the glc grid
    type(mct_aVect) :: glc_frac_this_ec_g

    ! attribute vector holding glc fraction in one elev class, on the land grid
    type(mct_aVect) :: glc_frac_this_ec_l

    ! attribute vector holding the product of (glc fraction in one elev class) x
    ! (icemask), on the glc grid
    type(mct_aVect) :: glc_frac_this_ec_times_icemask_g
    
    ! attribute vector holding fields to map (other than fraction) in one elevation
    ! class, on the land grid
    type(mct_aVect) :: glc_fields_this_ec_l
    
    character(len=*), parameter :: subname = 'map_glc2lnd_ec'
    !-----------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Determine attribute vector sizes
    ! ------------------------------------------------------------------------

    lsize_g = mct_aVect_lsize(g2x_g)
    lsize_l = mct_aVect_lsize(g2x_l)

    ! ------------------------------------------------------------------------
    ! Extract special fields from g2x_g
    ! ------------------------------------------------------------------------

    allocate(glc_frac(lsize_g))
    allocate(glc_topo(lsize_g))
    call mct_aVect_exportRattr(g2x_g, frac_field, glc_frac)
    call mct_aVect_exportRattr(g2x_g, topo_field, glc_topo)
    
    ! ------------------------------------------------------------------------
    ! Determine elevation class of each glc point
    ! ------------------------------------------------------------------------

    allocate(glc_elevclass(lsize_g))
    allocate(glc_frac_this_ec(lsize_g))
    call get_glc_elevation_classes(glc_topo, glc_elevclass)

    ! ------------------------------------------------------------------------
    ! Map each elevation class
    ! ------------------------------------------------------------------------

    call shr_string_listMerge(extra_fields, topo_field, fields_to_map)
    
    do n = 0, glc_get_num_elevation_classes()

       ! ------------------------------------------------------------------------
       ! Put fraction in this elevation class into an attribute vector
       ! ------------------------------------------------------------------------

       call get_frac_this_ec(glc_frac, glc_elevclass, n, glc_frac_this_ec)
       call mct_aVect_init(glc_frac_this_ec_g, rList = frac_field, lsize = lsize_g)
       call mct_aVect_importRattr(glc_frac_this_ec_g, frac_field, glc_frac_this_ec)

       ! ------------------------------------------------------------------------
       ! Map fraction to the land grid
       ! ------------------------------------------------------------------------

       call mct_aVect_init(glc_frac_this_ec_l, rList = frac_field, lsize = lsize_l)

       call seq_map_map(mapper = mapper, av_s = glc_frac_this_ec_g, av_d = glc_frac_this_ec_l, &
            norm = .true., avwts_s = g2x_g, avwtsfld_s = icemask_field)

       elevclass_as_string = glc_elevclass_as_string(n)
       frac_field_ec = frac_field // elevclass_as_string
       call mct_aVect_copy(glc_frac_this_ec_l, g2x_l, &
            rList = frac_field, TrList = frac_field_ec)
       
       ! ------------------------------------------------------------------------
       ! Map other fields to the land grid
       ! 
       ! Note that bare land values are mapped in the same way as ice-covered values
       ! ------------------------------------------------------------------------

       ! Create a mask that is (fraction in this elevation class) x (icemask). So, only
       ! grid cells that are both (a) within the icemask and (b) in this elevation class
       ! will be included in the following mapping.
       call make_aVect_frac_times_icemask(frac_av = glc_frac_this_ec_g, &
                                          mask_av = g2x_g, &
                                          frac_field = frac_field, &
                                          icemask_field = icemask_field, &
                                          frac_times_icemask_av = glc_frac_this_ec_times_icemask_g)

       call mct_aVect_init(glc_fields_this_ec_l, rList = fields_to_map, lsize = lsize_l)
       call seq_map_map(mapper = mapper, av_s = g2x_g, av_d = glc_fields_this_ec_l, &
            fldlist = fields_to_map, &
            norm = .true., &
            avwts_s = glc_frac_this_ec_times_icemask_g, &
            avwtsfld_s = frac_times_icemask_field)

       call set_topo_in_virtual_columns(n, glc_frac_this_ec_l, &
            frac_field, topo_field, &
            glc_fields_this_ec_l)

       call shr_string_listAddSuffix(fields_to_map, glc_elevclass_as_string(n), fields_to_map_ec)
       call mct_aVect_copy(glc_fields_this_ec_l, g2x_l, &
            rList = fields_to_map, TrList = fields_to_map_ec)

       ! ------------------------------------------------------------------------
       ! Clean up
       ! ------------------------------------------------------------------------
       
       call mct_aVect_clean(glc_frac_this_ec_l)
       call mct_aVect_clean(glc_frac_this_ec_g)
       call mct_aVect_clean(glc_frac_this_ec_times_icemask_g)
       call mct_aVect_clean(glc_fields_this_ec_l)

    end do

    deallocate(glc_frac)
    deallocate(glc_topo)
    deallocate(glc_frac_this_ec)
    
  end subroutine map_glc2lnd_ec


  !-----------------------------------------------------------------------
  subroutine get_glc_elevation_classes(glc_topo, glc_elevclass)
    !
    ! !DESCRIPTION:
    ! Get elevation class of each grid cell on the glc grid.
    !
    ! This does not consider glc_frac: it simply gives the elevation class that the grid
    ! cell would be in if it were ice-covered. So it never returns an elevation class of
    ! 0 (bare land). (This design would allow us, in the future, to have glc grid cells
    ! that are part ice-covered, part ice-free.)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: glc_topo(:)      ! topographic height
    integer , intent(out) :: glc_elevclass(:) ! elevation class
    !
    ! !LOCAL VARIABLES:
    integer :: npts
    integer :: glc_pt
    integer :: err_code
    
    character(len=*), parameter :: subname = 'get_glc_elevation_classes'
    !-----------------------------------------------------------------------

    npts = size(glc_elevclass)
    SHR_ASSERT((size(glc_topo) == npts), errMsg(__FILE__, __LINE__))

    do glc_pt = 1, npts
       call glc_get_elevation_class(glc_topo(glc_pt), glc_elevclass(glc_pt), err_code)
       select case (err_code)
       case (GLC_ELEVCLASS_ERR_NONE)
          ! Do nothing
       case (GLC_ELEVCLASS_ERR_TOO_LOW, GLC_ELEVCLASS_ERR_TOO_HIGH)
          write(logunit,*) subname, ': WARNING, for glc_pt, topo = ', glc_pt, glc_topo(glc_pt)
          write(logunit,*) glc_errcode_to_string(err_code)
       case default
          write(logunit,*) subname, ': ERROR getting elevation class for glc_pt = ', glc_pt
          write(logunit,*) glc_errcode_to_string(err_code)
          call shr_sys_abort(subname//': ERROR getting elevation class')
       end select
    end do
       
  end subroutine get_glc_elevation_classes

  !-----------------------------------------------------------------------
  subroutine get_frac_this_ec(glc_frac, glc_elevclass, this_elevclass, glc_frac_this_ec)
    !
    ! !DESCRIPTION:
    ! Get fractional ice coverage in a given elevation class.
    !
    ! The input glc_elevclass gives the elevation class of each glc grid cell, assuming
    ! that the grid cell is ice-covered.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: glc_frac(:)         ! total ice sheet fraction in each glc grid cell
    integer , intent(in)  :: glc_elevclass(:)    ! elevation class of each glc grid cell
    integer , intent(in)  :: this_elevclass      ! elevation class index of interest
    real(r8), intent(out) :: glc_frac_this_ec(:) ! ice fraction in this elevation class
    !
    ! !LOCAL VARIABLES:
    integer :: npts
    
    character(len=*), parameter :: subname = 'get_frac_this_ec'
    !-----------------------------------------------------------------------

    npts = size(glc_frac_this_ec)
    SHR_ASSERT((size(glc_frac) == npts), errMsg(__FILE__, __LINE__))
    SHR_ASSERT((size(glc_elevclass) == npts), errMsg(__FILE__, __LINE__))
    
    if (this_elevclass == 0) then
       glc_frac_this_ec(:) = 1._r8 - glc_frac(:)
    else
       where (glc_elevclass == this_elevclass)
          glc_frac_this_ec = glc_frac
       elsewhere
          glc_frac_this_ec = 0._r8
       end where
    end if
       
  end subroutine get_frac_this_ec

  !-----------------------------------------------------------------------
  subroutine set_topo_in_virtual_columns(elev_class, glc_frac_this_ec_l, &
       frac_field, topo_field, &
       glc_topo_this_ec_l)
    !
    ! !DESCRIPTION:
    ! Sets the topo field for virtual columns, in a given elevation class.
    !
    ! This is needed because virtual columns (i.e., elevation classes that have no
    ! contributing glc grid cells) won't have any topographic information mapped onto
    ! them, so would otherwise end up with an elevation of 0.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: elev_class
    type(mct_aVect), intent(in) :: glc_frac_this_ec_l  ! attr vect containing frac_field
    character(len=*), intent(in) :: frac_field
    character(len=*), intent(in) :: topo_field
    type(mct_aVect), intent(inout) :: glc_topo_this_ec_l  ! attr vect containing topo_field
    !
    ! !LOCAL VARIABLES:
    integer :: lsize
    real(r8) :: topo_virtual

    ! The following need to be pointers to satisfy the MCT interface:
    real(r8), pointer :: frac_l(:)  ! ice fraction in this elev class, land grid
    real(r8), pointer :: topo_l(:)  ! topographic height in this elev class, land grid
    
    character(len=*), parameter :: subname = 'set_virtual_elevation_classes'
    !-----------------------------------------------------------------------

    ! Extract fields from attribute vectors
    lsize = mct_aVect_lsize(glc_frac_this_ec_l)
    SHR_ASSERT(mct_aVect_lsize(glc_topo_this_ec_l) == lsize, errMsg(__FILE__, __LINE__))
    allocate(frac_l(lsize))
    allocate(topo_l(lsize))
    call mct_aVect_exportRattr(glc_frac_this_ec_l, frac_field, frac_l)
    call mct_aVect_exportRattr(glc_topo_this_ec_l, topo_field, topo_l)

    ! Set topo field for virtual columns
    topo_virtual = glc_mean_elevation_virtual(elev_class)
    where (frac_l <= 0)
       topo_l = topo_virtual
    end where
    
    ! Put updated field back in attribute vector
    call mct_aVect_importRattr(glc_topo_this_ec_l, topo_field, topo_l)
    
    deallocate(frac_l)
    deallocate(topo_l)
    
  end subroutine set_topo_in_virtual_columns

  !-----------------------------------------------------------------------
  subroutine make_aVect_frac_times_icemask(frac_av, mask_av, frac_field, icemask_field, &
       frac_times_icemask_av)
    !
    ! !DESCRIPTION:
    ! Create an attribute vector that is the product of frac_field and icemask_field
    !
    ! The resulting frac_times_icemask_av will have a field frac_times_icemask_field which
    ! contains this product. This attribute vector is initialized here; it is expected to
    ! come in in an uninitialized/cleaned state. (So it needs to be cleaned with a call to
    ! mct_aVect_clean later - including before the next call to this routine.)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect), intent(in)  :: frac_av  ! attr vect containing frac_field
    type(mct_aVect), intent(in)  :: mask_av  ! attr vect containing icemask_field
    character(len=*), intent(in) :: frac_field
    character(len=*), intent(in) :: icemask_field
    type(mct_aVect), intent(out) :: frac_times_icemask_av  ! attr vect that will contain frac_times_icemask_field
    !
    ! !LOCAL VARIABLES:
    integer :: lsize
    
    character(len=*), parameter :: subname = 'make_aVect_frac_times_icemask'
    !-----------------------------------------------------------------------

    lsize = mct_aVect_lsize(frac_av)
    SHR_ASSERT(mct_aVect_lsize(mask_av) == lsize, errMsg(__FILE__, __LINE__))

    call mct_aVect_init(frac_times_icemask_av, rList = frac_times_icemask_field, lsize = lsize)
    call mct_aVect_copy(aVin = frac_av, aVout = frac_times_icemask_av, &
         rList = frac_field, TrList = frac_times_icemask_field)
    call mct_aVect_mult(frac_times_icemask_av, mask_av, icemask_field)
    
  end subroutine make_aVect_frac_times_icemask
  
end module map_glc2lnd_mod

  
