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

  public :: map_glc2lnd_ec ! moab (tag-based) version of map_glc2lnd_ec
  public :: map_glc2lnd_ec_init ! define the moab tags map_glc2lnd_ec needs

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: get_glc_elevation_classes     ! get elevation class of each glc cell
  private :: get_frac_this_ec              ! get fraction in a given elevation class

  ! scratch tag holding M(icemask), the denominator of the per-EC fraction
  ! normalization in map_glc2lnd_ec; needed on both the glc and lnd meshes
  character(len=*), parameter :: scratch_icemask_num = 'Sg_icemsk_num'

contains

  !-----------------------------------------------------------------------
  subroutine map_glc2lnd_ec(mapper, &
       frac_field, topo_field, icemask_field, extra_fields, frac_l_out)
    !
    ! !DESCRIPTION:
    ! Maps fields from the GLC mesh (mapper%src_mbid)
    ! to the LND mesh (mapper%tgt_mbid), separated by elevation class, operating on
    ! moab tags. The per-EC results are written into the tags named
    ! <field><EC> on the land mesh (these are members of seq_flds_x2l_fields, already
    ! defined there).
    !
    ! The driver-mct version does, for each elevation class n, two normalized maps:
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

    character(len=*), parameter :: subname = 'map_glc2lnd_ec'
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
          ! The divisions are done as reciprocal multiplies (with a zero denominator
          ! passing through as a multiply by zero), exactly like the normalization in
          ! the mct seq_map_avNormArr, so the results match the mct driver bit for bit.
          ! frac_n_l = M(w_n) / M(icemask)
          denom = num_l(i, kmask)
          if (denom /= 0.0_r8) then
             denom = 1.0_r8/denom
          end if
          out_l(i, kfrac0) = num_l(i, kfrac0) * denom
          ! field_n_l = M(field*w_n) / M(w_n)
          denom = num_l(i, kfrac0)
          if (denom /= 0.0_r8) then
             denom = 1.0_r8/denom
          end if
          out_l(i, ktopo0) = num_l(i, ktopo0) * denom
          do k = 1, nextra
             out_l(i, ktopo0+k) = num_l(i, ktopo0+k) * denom
          end do
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

  end subroutine map_glc2lnd_ec

  !-----------------------------------------------------------------------
  subroutine map_glc2lnd_ec_init(glc_mbid, lnd_mbid)
    !
    ! !DESCRIPTION:
    ! Defines the moab tags that map_glc2lnd_ec needs, beyond the per-EC
    ! destination tags that already exist on the land mesh as part of the x2l
    ! field list: the same per-EC names on the glc mesh, where they stage the
    ! pre-multiplied numerators, and the scratch tag for the icemask
    ! normalization, which is needed on both meshes.
    !
    ! Called from whichever prep module sets up the glc->lnd map: prep_lnd when
    ! glc feeds back to land, prep_glc otherwise (the smb renormalization goes
    ! through the same mapping).
    !
    ! !USES:
    use iMOAB, only : iMOAB_DefineTagStorage
    use iso_c_binding, only : C_NULL_CHAR
    use seq_flds_mod, only : seq_flds_x2l_fields_from_glc
    !
    ! !ARGUMENTS:
    integer, intent(in) :: glc_mbid  ! iMOAB app id of the glc mesh on the coupler
    integer, intent(in) :: lnd_mbid  ! iMOAB app id of the lnd mesh on the coupler
    !
    ! !LOCAL VARIABLES:
    integer :: ierr, tagtype, numco, tagindex
    character(CXX) :: tagname

    character(len=*), parameter :: subname = 'map_glc2lnd_ec_init'
    !-----------------------------------------------------------------------

    if (glc_mbid < 0 .or. lnd_mbid < 0) return

    tagtype = 1  ! dense, double
    numco = 1

    tagname = trim(seq_flds_x2l_fields_from_glc)//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(glc_mbid, tagname, tagtype, numco, tagindex)
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR defining glc per-EC staging tags on the glc mesh')
    endif

    tagname = scratch_icemask_num//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(glc_mbid, tagname, tagtype, numco, tagindex)
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR defining '//scratch_icemask_num//' on the glc mesh')
    endif
    ierr = iMOAB_DefineTagStorage(lnd_mbid, tagname, tagtype, numco, tagindex)
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR defining '//scratch_icemask_num//' on the land mesh')
    endif

  end subroutine map_glc2lnd_ec_init

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

end module map_glc2lnd_mod
