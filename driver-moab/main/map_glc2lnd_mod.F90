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
  use shr_kind_mod, only : cxx => SHR_KIND_CXX
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
  public :: map_glc2lnd_ec_moab ! moab (tag-based) version of map_glc2lnd_ec

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
  subroutine map_glc2lnd_ec_moab(mapper, &
       frac_field, topo_field, icemask_field, extra_fields, frac_l_out)
    !
    ! !DESCRIPTION:
    ! moab version of map_glc2lnd_ec: maps fields from the GLC mesh (mapper%src_mbid)
    ! to the LND mesh (mapper%tgt_mbid), separated by elevation class, operating on
    ! moab tags. The per-EC results are written into the tags named
    ! <field><EC> on the land mesh (these are members of seq_flds_x2l_fields, already
    ! defined there).
    !
    ! The mct version does, for each elevation class n, two normalized maps:
    !     frac_n_l  = M(icemask*frac_n) / M(icemask)
    !     field_n_l = M(icemask*frac_n*field) / M(icemask*frac_n)
    ! Because the same conservative map M is used throughout, this is reproduced here
    ! with ONE raw (norm=.false.) map of all pre-multiplied numerators, followed by
    ! explicit divisions on the land side (the raw row sums cancel in the ratios).
    ! The numerators are staged in the per-EC destination tag names on the glc mesh
    ! (defined at init), plus one scratch tag 'Sg_icemsk_num' for M(icemask).
    !
    ! Assumes the following tags exist:
    ! - on the glc mesh: frac_field, topo_field, icemask_field, each extra field,
    !   the per-EC names <field><EC> for frac/topo/extras, and Sg_icemsk_num
    ! - on the lnd mesh: the same per-EC names and Sg_icemsk_num
    !
    ! !USES:
    use iMOAB, only : iMOAB_GetMeshInfo, iMOAB_GetDoubleTagStorage, iMOAB_SetDoubleTagStorage
    use iso_c_binding, only : C_NULL_CHAR
    use glc_elevclass_mod, only : glc_all_elevclass_strings, GLC_ELEVCLASS_STRLEN
    !
    ! !ARGUMENTS:
    type(seq_map), intent(inout) :: mapper        ! conservative glc->lnd mapper with moab context
    character(len=*), intent(in) :: frac_field    ! name of glc field containing glc ice fraction
    character(len=*), intent(in) :: topo_field    ! name of glc field containing glc topo
    character(len=*), intent(in) :: icemask_field ! name of glc field containing ice mask
    character(len=*), intent(in) :: extra_fields  ! colon-delimited additional fields ('' or ' ' for none)
    real(r8), optional, intent(out) :: frac_l_out(:,0:) ! normalized per-EC fractions on the land mesh
    !
    ! !LOCAL VARIABLES:
    integer :: lsize_g, lsize_l
    integer :: nEC, nextra, nnum
    integer :: n, i, k, kfrac0, ktopo0, kextra0, kmask
    integer :: ierr, ent_type, arrsize
    integer :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)

    real(r8), allocatable :: glc_frac(:)     ! total ice fraction in each glc cell
    real(r8), allocatable :: glc_topo(:)     ! topographic height of each glc cell
    real(r8), allocatable :: glc_icemask(:)  ! icemask of each glc cell
    real(r8), allocatable :: glc_extra(:,:)  ! extra fields on the glc mesh
    real(r8), allocatable :: frac_this_ec(:) ! ice fraction in one elevation class
    integer , allocatable :: glc_elevclass(:)
    real(r8), allocatable :: num_g(:,:)      ! numerators on the glc mesh
    real(r8), allocatable :: num_l(:,:)      ! mapped numerators on the land mesh
    real(r8), allocatable :: out_l(:,:)      ! normalized per-EC fields on the land mesh
    real(r8) :: denom, topo_virtual

    character(len=GLC_ELEVCLASS_STRLEN), allocatable :: ec_strings(:)
    character(CXX) :: numlist   ! colon-separated list of all numerator tag names
    character(CXX) :: tagname
    type(mct_list) :: temp_list
    type(mct_string) :: mctOStr
    character(CXX) :: extra_names(20)  ! names of the extra fields (assumed few)

    ! dummy zero-size attribute vectors to satisfy the seq_map_map interface
    type(mct_aVect) :: av_dum_s, av_dum_d

    character(len=*), parameter :: scratch_icemask_num = 'Sg_icemsk_num'
    character(len=*), parameter :: subname = 'map_glc2lnd_ec_moab'
    !-----------------------------------------------------------------------

    if (mapper%src_mbid < 0 .or. mapper%tgt_mbid < 0) return

    ent_type = 1 ! cells on both meshes

    ierr  = iMOAB_GetMeshInfo ( mapper%src_mbid, nvert, nvise, nbl, nsurf, nvisBC )
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting glc mesh info')
    lsize_g = nvise(1)
    ierr  = iMOAB_GetMeshInfo ( mapper%tgt_mbid, nvert, nvise, nbl, nsurf, nvisBC )
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting lnd mesh info')
    lsize_l = nvise(1)

    nEC = glc_get_num_elevation_classes()

    ! parse the extra field names
    nextra = 0
    if (len_trim(extra_fields) > 0) then
       call mct_list_init(temp_list, extra_fields)
       nextra = mct_list_nitem(temp_list)
       if (nextra > size(extra_names)) call shr_sys_abort(subname//' ERROR too many extra fields')
       do i = 1, nextra
          call mct_list_get(mctOStr, i, temp_list)
          extra_names(i) = mct_string_toChar(mctOStr)
          call mct_string_clean(mctOStr)
       end do
       call mct_list_clean(temp_list)
    end if

    ! ------------------------------------------------------------------------
    ! Extract needed fields from the glc mesh tags
    ! ------------------------------------------------------------------------

    allocate(glc_frac(lsize_g), glc_topo(lsize_g), glc_icemask(lsize_g))
    allocate(frac_this_ec(lsize_g), glc_elevclass(lsize_g))
    tagname = trim(frac_field)//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mapper%src_mbid, tagname, lsize_g, ent_type, glc_frac)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting '//trim(frac_field))
    tagname = trim(topo_field)//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mapper%src_mbid, tagname, lsize_g, ent_type, glc_topo)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting '//trim(topo_field))
    tagname = trim(icemask_field)//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mapper%src_mbid, tagname, lsize_g, ent_type, glc_icemask)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting '//trim(icemask_field))
    if (nextra > 0) then
       allocate(glc_extra(lsize_g, nextra))
       do i = 1, nextra
          tagname = trim(extra_names(i))//C_NULL_CHAR
          ierr = iMOAB_GetDoubleTagStorage(mapper%src_mbid, tagname, lsize_g, ent_type, glc_extra(:,i))
          if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting '//trim(extra_names(i)))
       end do
    end if

    call get_glc_elevation_classes(glc_topo, glc_elevclass)

    ! ------------------------------------------------------------------------
    ! Build the numerator fields and the combined tag list
    ! layout per elevation class n (0..nEC):
    !   <frac_field>NN  = frac_n * icemask                      (this is M's input w_n)
    !   <topo_field>NN  = topo * w_n
    !   <extra>NN       = extra * w_n
    ! plus one extra column: Sg_icemsk_num = icemask
    ! ------------------------------------------------------------------------

    nnum = (nEC+1)*(2+nextra) + 1
    allocate(num_g(lsize_g, nnum))
    allocate(ec_strings(0:nEC))
    ec_strings = glc_all_elevclass_strings(include_zero = .true.)

    numlist = ''
    k = 0
    kfrac0 = 1  ! columns kfrac0 + n*(2+nextra) hold w_n, etc.
    do n = 0, nEC
       call get_frac_this_ec(glc_frac, glc_elevclass, n, frac_this_ec)
       ! frac numerator: w_n = frac_n * icemask
       k = k + 1
       num_g(:,k) = frac_this_ec(:) * glc_icemask(:)
       call add_to_list(numlist, trim(frac_field)//trim(ec_strings(n)))
       ! topo numerator: topo * w_n
       k = k + 1
       num_g(:,k) = glc_topo(:) * num_g(:,k-1)
       call add_to_list(numlist, trim(topo_field)//trim(ec_strings(n)))
       ! extra numerators
       do i = 1, nextra
          k = k + 1
          num_g(:,k) = glc_extra(:,i) * num_g(:,k-1-i)
          call add_to_list(numlist, trim(extra_names(i))//trim(ec_strings(n)))
       end do
    end do
    ! icemask numerator (for the frac normalization)
    k = k + 1
    kmask = k
    num_g(:,k) = glc_icemask(:)
    call add_to_list(numlist, scratch_icemask_num)

    ! set the numerators on the glc mesh tags
    arrsize = nnum * lsize_g
    tagname = trim(numlist)//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage(mapper%src_mbid, tagname, arrsize, ent_type, num_g)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR setting numerator tags on glc mesh')

    ! ------------------------------------------------------------------------
    ! One raw map of all numerators glc -> lnd
    ! ------------------------------------------------------------------------

    call mct_aVect_init(av_dum_s, rList = frac_field, lsize = 0)
    call mct_aVect_init(av_dum_d, rList = frac_field, lsize = 0)
    call seq_map_map(mapper, av_dum_s, av_dum_d, fldlist=trim(numlist), norm=.false.)
    call mct_aVect_clean(av_dum_s)
    call mct_aVect_clean(av_dum_d)

    ! ------------------------------------------------------------------------
    ! Get mapped numerators on the land mesh and normalize
    ! ------------------------------------------------------------------------

    allocate(num_l(lsize_l, nnum))
    allocate(out_l(lsize_l, nnum-1))
    num_l = 0.0_r8
    arrsize = nnum * lsize_l
    ierr = iMOAB_GetDoubleTagStorage(mapper%tgt_mbid, tagname, arrsize, ent_type, num_l)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting numerator tags on lnd mesh')

    do n = 0, nEC
       kfrac0 = 1 + n*(2+nextra)   ! column of w_n
       ktopo0 = kfrac0 + 1
       topo_virtual = glc_mean_elevation_virtual(n)
       do i = 1, lsize_l
          ! frac_n_l = M(w_n) / M(icemask)
          denom = num_l(i, kmask)
          if (denom /= 0.0_r8) then
             out_l(i, kfrac0) = num_l(i, kfrac0) / denom
          else
             out_l(i, kfrac0) = 0.0_r8
          end if
          ! field_n_l = M(field*w_n) / M(w_n)
          denom = num_l(i, kfrac0)
          if (denom /= 0.0_r8) then
             out_l(i, ktopo0) = num_l(i, ktopo0) / denom
             do k = 1, nextra
                out_l(i, ktopo0+k) = num_l(i, ktopo0+k) / denom
             end do
          else
             out_l(i, ktopo0) = 0.0_r8
             do k = 1, nextra
                out_l(i, ktopo0+k) = 0.0_r8
             end do
          end if
          ! set the topo field for virtual columns (no contributing glc cells)
          if (out_l(i, kfrac0) <= 0.0_r8) then
             out_l(i, ktopo0) = topo_virtual
          end if
       end do
       if (present(frac_l_out)) then
          frac_l_out(:, n) = out_l(:, kfrac0)
       end if
    end do

    ! ------------------------------------------------------------------------
    ! Store the normalized per-EC fields back into the land mesh tags
    ! (all list entries except the trailing Sg_icemsk_num scratch)
    ! ------------------------------------------------------------------------

    i = len_trim(numlist) - len(scratch_icemask_num) - 1 ! strip ':Sg_icemsk_num'
    tagname = numlist(1:i)//C_NULL_CHAR
    arrsize = (nnum-1) * lsize_l
    ierr = iMOAB_SetDoubleTagStorage(mapper%tgt_mbid, tagname, arrsize, ent_type, out_l)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR setting per-EC tags on lnd mesh')

    deallocate(glc_frac, glc_topo, glc_icemask, frac_this_ec, glc_elevclass)
    if (allocated(glc_extra)) deallocate(glc_extra)
    deallocate(num_g, num_l, out_l, ec_strings)

  contains

    subroutine add_to_list(list, name)
      ! append name to a colon-separated list
      character(len=*), intent(inout) :: list
      character(len=*), intent(in)    :: name
      if (len_trim(list) == 0) then
         list = trim(name)
      else
         list = trim(list)//':'//trim(name)
      end if
    end subroutine add_to_list

  end subroutine map_glc2lnd_ec_moab

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
    SHR_ASSERT_FL((size(glc_topo) == npts), __FILE__, __LINE__)

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
    SHR_ASSERT_FL((size(glc_frac) == npts), __FILE__, __LINE__)
    SHR_ASSERT_FL((size(glc_elevclass) == npts), __FILE__, __LINE__)

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
    SHR_ASSERT_FL(mct_aVect_lsize(glc_topo_this_ec_l) == lsize, __FILE__, __LINE__)
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
    SHR_ASSERT_FL(mct_aVect_lsize(mask_av) == lsize, __FILE__, __LINE__)

    call mct_aVect_init(frac_times_icemask_av, rList = frac_times_icemask_field, lsize = lsize)
    call mct_aVect_copy(aVin = frac_av, aVout = frac_times_icemask_av, &
         rList = frac_field, TrList = frac_times_icemask_field)
    call mct_aVect_mult(frac_times_icemask_av, mask_av, icemask_field)

  end subroutine make_aVect_frac_times_icemask

end module map_glc2lnd_mod
