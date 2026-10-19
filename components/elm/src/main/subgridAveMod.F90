module subgridAveMod

#include "shr_assert.h"
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Utilities to perfrom subgrid averaging
  !
  ! Each average is computed by one OpenACC loop over the parent level (an
  ! independent weighted sum per parent), with the children summed in
  ! ascending index order (loop seq), which is the order of the serial
  ! scatter loops these replace, so the results are bit-for-bit with them.
  ! Children of a landunit or column are contiguous (lun_pp%coli/colf,
  ! lun_pp%pfti/pftf, col_pp%pfti/pftf); children of a topounit or gridcell
  ! are not (see elm_ptrs_compdown), so those use the child maps built in
  ! init_subgrid_ave.
  !
  ! Averages to gridcells (p2g, c2g, l2g, t2g) call the averages of
  ! subgridAveViewMod (one vector routine per chunk of gridcells) with
  ! elm_subgrid_view, the subgrid structure set up in init_subgrid_ave; the
  ! history uses the same averages.
  !
  ! The scale factors come from lookup tables built once in init_subgrid_ave:
  ! main_scale_c2l(c, c2l_scale_type) and main_scale_l2g_lookup(l2g_scale_type,
  ! lun_pp%itype(l)) (also used for l2t). The averaging routines take the
  ! integer scale types below; the generic interfaces also accept the scale
  ! type names ('unity', 'urbanf', ...).
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use column_varcon , only : icol_roof, icol_sunwall, icol_shadewall
  use column_varcon , only : icol_road_perv , icol_road_imperv
  use elm_varcon    , only : grlnd, nameg, namet, namel, namec, namep,spval
  use elm_varctl    , only : iulog
  use decompMod     , only : bounds_type
  use GridcellType  , only : grc_pp
  use TopounitType  , only : top_pp
  use LandunitType  , only : lun_pp
  use ColumnType    , only : col_pp
  use VegetationType, only : veg_pp
  use abortutils    , only : endrun
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  use subgridAveViewMod, only : subgrid_ave_view, view_spval => spval
  use subgridAveViewMod, only : p2g_range_unity, c2g_range_unity
  use subgridAveViewMod, only : build_child_map, subgrid_ave_view_attach, subgrid_ave_view_detach
  use subgridAveViewMod, only : p2g_range, c2g_range, l2g_range, t2g_range
  !
  ! !PUBLIC TYPES:
  implicit none
  save

  ! l2g (and l2t) scale types: index into main_scale_l2g_lookup
  integer, public, parameter :: unity = 1, natveg = 2, veg = 3
  integer, public, parameter :: ice = 4, nonurb = 5, lake = 6
  ! c2l scale types (unity = 1 as well): index into main_scale_c2l
  integer, public, parameter :: urbanf = 2, urbans = 3
  integer, parameter :: num_c2l_types = 3, num_l2g_types = 6
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_subgrid_ave ! Build the scale lookups and child maps
  public :: p2c   ! Perform an average pfts to columns
  public :: p2l   ! Perform an average pfts to landunits
  public :: p2g   ! Perform an average pfts to gridcells
  public :: c2l   ! Perform an average columns to landunits
  public :: c2g   ! Perform an average columns to gridcells
  public :: l2g   ! Perform an average landunits to gridcells
  public :: p2t   ! Perform an average pfts to topounits
  public :: c2t   ! Perform an average columns to topounits
  public :: l2t   ! Perform an average landunits to topounits
  public :: t2g   ! Perform an average topounits to gridcells
  public :: p2c_1d, p2c_1d_filter, p2l_1d

  ! Subgrid structure for averages to gridcells (see subgridAveViewMod)
  type(subgrid_ave_view), public, protected :: elm_subgrid_view
  logical, private :: view_on_device = .false.
  ! gridcells per gang in the averages to gridcells (each gang runs one
  ! vector routine over its chunk)
  integer, parameter, private :: gridcell_chunk = 128

  real(r8), parameter, private :: sumwt_tol = 1.0_r8 + 1.e-6_r8 ! sumwt above this is an error

  interface p2c
     module procedure p2c_1d
     module procedure p2c_1d_name
     module procedure p2c_2d
     module procedure p2c_2d_name
     module procedure p2c_1d_filter
  end interface
  interface p2l
     module procedure p2l_1d
     module procedure p2l_1d_name
     module procedure p2l_2d
     module procedure p2l_2d_name
  end interface
  interface p2g
     module procedure p2g_1d
     module procedure p2g_1d_name
     module procedure p2g_2d
     module procedure p2g_2d_name
  end interface
  interface c2l
     module procedure c2l_1d
     module procedure c2l_1d_name
     module procedure c2l_2d
     module procedure c2l_2d_name
  end interface
  interface c2g
     module procedure c2g_1d
     module procedure c2g_1d_name
     module procedure c2g_2d
     module procedure c2g_2d_name
  end interface
  interface l2g
     module procedure l2g_1d
     module procedure l2g_1d_name
     module procedure l2g_2d
     module procedure l2g_2d_name
  end interface
  interface p2t
     module procedure p2t_1d
     module procedure p2t_1d_name
     module procedure p2t_2d
     module procedure p2t_2d_name
  end interface
  interface c2t
     module procedure c2t_1d
     module procedure c2t_1d_name
     module procedure c2t_2d
     module procedure c2t_2d_name
  end interface
  interface l2t
     module procedure l2t_1d
     module procedure l2t_1d_name
     module procedure l2t_2d
     module procedure l2t_2d_name
  end interface
  interface t2g
     module procedure t2g_1d
     module procedure t2g_2d
  end interface

  ! Scale factor lookups (histFileMod copies them into the history scale tables)
  real(r8), allocatable, public, protected, target :: main_scale_l2g_lookup(:,:) ! (l2g_scale_type, lun_pp%itype)
  !$acc declare create(main_scale_l2g_lookup(:,:))
  real(r8), allocatable, public, protected, target :: main_scale_c2l(:,:)        ! (column, c2l_scale_type)
  !$acc declare create(main_scale_c2l(:,:))

  ! Child maps: the children of parent k are idx(ptr(k):ptr(k+1)-1), in ascending order
  ! (landunits of a gridcell, for elm_subgrid_view; points of a topounit, for x2t)
  integer, allocatable, private, target :: grc_lun_ptr(:), grc_lun_idx(:)
  integer, allocatable, private :: top_pft_ptr(:), top_pft_idx(:)
  integer, allocatable, private :: top_col_ptr(:), top_col_idx(:)
  integer, allocatable, private :: top_lun_ptr(:), top_lun_idx(:)
  !$acc declare create(grc_lun_ptr(:), grc_lun_idx(:), top_pft_ptr(:), top_pft_idx(:), &
  !$acc                top_col_ptr(:), top_col_idx(:), top_lun_ptr(:), top_lun_idx(:))

  ! WJS (10-14-11): TODO:
  !
  ! - I believe that scale_p2c, scale_c2l and scale_l2g should be included in the sumwt
  ! accumulations (e.g., sumwt = sumwt + wtgcell * scale_p2c * scale_c2l * scale_l2g), but
  ! that requires some more thought to (1) make sure that is correct, and (2) make sure it
  ! doesn't break the urban scaling. (See also my notes in init_scale_l2g_lookup.)
  !   - Once that is done, you could use a scale of 0, avoiding the need for the use of
  !   spval and the special checks that requires.
  ! -----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine init_subgrid_ave(bounds)
    ! !DESCRIPTION:
    ! Build the scale factor lookups, the topounit/gridcell child maps and
    ! elm_subgrid_view. Call once with the processor bounds, after the subgrid
    ! structure and lun_pp%canyon_hwr are set (call again if canyon_hwr changes).
    ! With OpenACC, the subgrid types' arrays must be on the device first.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds ! processor bounds
    !-----------------------------------------------------------------------

    if (view_spval /= spval) then
       write(iulog,*) 'init_subgrid_ave error: subgridAveViewMod spval differs from elm_varcon spval'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    if (view_on_device) call subgrid_ave_view_detach(elm_subgrid_view)

    call init_scale_l2g_lookup()
    call init_scale_c2l(bounds)

    associate(begg => bounds%begg, endg => bounds%endg, begt => bounds%begt, endt => bounds%endt, &
              begl => bounds%begl, endl => bounds%endl, begc => bounds%begc, endc => bounds%endc, &
              begp => bounds%begp, endp => bounds%endp)

    if (allocated(grc_lun_ptr)) then
       deallocate(grc_lun_ptr, grc_lun_idx)
       deallocate(top_pft_ptr, top_pft_idx, top_col_ptr, top_col_idx, top_lun_ptr, top_lun_idx)
    end if
    allocate(grc_lun_ptr(begg:endg+1), grc_lun_idx(begl:endl))
    allocate(top_pft_ptr(begt:endt+1), top_pft_idx(begp:endp))
    allocate(top_col_ptr(begt:endt+1), top_col_idx(begc:endc))
    allocate(top_lun_ptr(begt:endt+1), top_lun_idx(begl:endl))

    call child_map(begg, endg, begl, endl, lun_pp%gridcell(begl:endl), grc_lun_ptr, grc_lun_idx)
    call child_map(begt, endt, begp, endp, veg_pp%topounit(begp:endp), top_pft_ptr, top_pft_idx)
    call child_map(begt, endt, begc, endc, col_pp%topounit(begc:endc), top_col_ptr, top_col_idx)
    call child_map(begt, endt, begl, endl, lun_pp%topounit(begl:endl), top_lun_ptr, top_lun_idx)

    end associate

    !$acc update device(grc_lun_ptr, grc_lun_idx, &
    !$acc               top_pft_ptr, top_pft_idx, top_col_ptr, top_col_idx, top_lun_ptr, top_lun_idx)

    call set_elm_subgrid_view(bounds)
    call subgrid_ave_view_attach(elm_subgrid_view)
    view_on_device = .true.

  end subroutine init_subgrid_ave

  !-----------------------------------------------------------------------
  subroutine child_map(begpar, endpar, begch, endch, parent, ptr, idx)
    ! !DESCRIPTION:
    ! build_child_map (subgridAveViewMod), aborting on a parent out of range.
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: begpar, endpar, begch, endch
    integer, intent(in)  :: parent(begch:endch)  ! parent index of each child
    integer, intent(out) :: ptr(begpar:endpar+1) ! children of k: idx(ptr(k):ptr(k+1)-1)
    integer, intent(out) :: idx(begch:endch)
    !
    ! !LOCAL VARIABLES:
    integer :: bad_child
    !-----------------------------------------------------------------------

    call build_child_map(begpar, endpar, begch, endch, parent, ptr, idx, bad_child)
    if (bad_child /= 0) then
       write(iulog,*) 'child_map error: parent ', parent(bad_child), ' of child ', bad_child, &
            ' outside ', begpar, endpar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine child_map

  !-----------------------------------------------------------------------
  subroutine set_elm_subgrid_view(bounds)
    ! !DESCRIPTION:
    ! Point elm_subgrid_view at the landunit index, the scale tables and the
    ! subgrid types' arrays (weights, active flags and snow layers are read
    ! through it at every average, so changes to them need no refresh).
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds ! processor bounds
    !-----------------------------------------------------------------------

    associate(v => elm_subgrid_view)
    v%begg = bounds%begg ; v%endg = bounds%endg
    v%begt = bounds%begt ; v%endt = bounds%endt
    v%begl = bounds%begl ; v%endl = bounds%endl
    v%begc = bounds%begc ; v%endc = bounds%endc
    v%begp = bounds%begp ; v%endp = bounds%endp
    v%grc_lun_ptr => grc_lun_ptr ; v%grc_lun_idx => grc_lun_idx
    v%topi => grc_pp%topi ; v%topf => grc_pp%topf
    v%coli => lun_pp%coli ; v%colf => lun_pp%colf
    v%pfti => lun_pp%pfti ; v%pftf => lun_pp%pftf
    v%wtgcell_t => top_pp%wtgcell ; v%active_t => top_pp%active
    v%wtgcell_l => lun_pp%wtgcell ; v%active_l => lun_pp%active
    v%wtgcell_c => col_pp%wtgcell ; v%active_c => col_pp%active
    v%wtgcell_p => veg_pp%wtgcell ; v%active_p => veg_pp%active
    v%col_of_pft => veg_pp%column
    v%lun_type   => lun_pp%itype
    v%scale_c2l        => main_scale_c2l
    v%scale_l2g_lookup => main_scale_l2g_lookup
    v%snl        => col_pp%snl
    end associate

  end subroutine set_elm_subgrid_view

  !-----------------------------------------------------------------------
  subroutine init_scale_l2g_lookup()
    ! DESCRIPTION:
    ! Build main_scale_l2g_lookup(l2g_scale_type, landunit type), the scale factor
    ! for each landunit type (used for both l2g and l2t averages).
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop, istice, istice_mec, istdlak
    use landunit_varcon, only : isturb_MIN, isturb_MAX, max_lunit
    !-----------------------------------------------------------------------

     ! ------------ WJS (10-14-11): IMPORTANT GENERAL NOTES ------------
     !
     ! Since scale_l2g is not currently included in the sumwt accumulations, you need to
     ! be careful about the scale values you use. Values of 1 and spval are safe
     ! (including having multiple landunits with value 1), but only use other values if
     ! you know what you are doing! For example, using a value of 0 is NOT the correct way
     ! to exclude a landunit from the average, because the normalization will be done
     ! incorrectly in this case: instead, use spval to exclude a landunit from the
     ! average. Similarly, using a value of 2 is NOT the correct way to give a landunit
     ! double relative weight in general, because the normalization won't be done
     ! correctly in this case, either.
     !
     ! In the longer-term, I believe that the correct solution to this problem is to
     ! include scale_l2g (and the other scale factors) in the sumwt accumulations
     ! (e.g., sumwt = sumwt + wtgcell * scale_p2c * scale_c2l * scale_l2g), but that
     ! requires some more thought to (1) make sure that is correct, and (2) make sure it
     ! doesn't break the urban scaling.
     !
     ! -----------------------------------------------------------------

     if (allocated(main_scale_l2g_lookup)) deallocate(main_scale_l2g_lookup)
     allocate(main_scale_l2g_lookup(num_l2g_types, max_lunit))

     ! Landunits that keep the default spval are excluded from the average
     main_scale_l2g_lookup(:,:) = spval
     main_scale_l2g_lookup(unity,:) = 1.0_r8
     main_scale_l2g_lookup(natveg,istsoil) = 1.0_r8
     main_scale_l2g_lookup(veg,istsoil) = 1.0_r8
     main_scale_l2g_lookup(veg,istcrop) = 1.0_r8
     main_scale_l2g_lookup(ice,istice) = 1.0_r8
     main_scale_l2g_lookup(ice,istice_mec) = 1.0_r8
     main_scale_l2g_lookup(nonurb,:) = 1.0_r8
     main_scale_l2g_lookup(nonurb,isturb_MIN:isturb_MAX) = spval
     main_scale_l2g_lookup(lake,istdlak) = 1.0_r8

     !$acc update device(main_scale_l2g_lookup(:,:))

  end subroutine init_scale_l2g_lookup

  !-----------------------------------------------------------------------
  subroutine init_scale_c2l(bounds)
    ! Description:
    ! Build main_scale_c2l(c, c2l_scale_type), the scale factor for each column
    ! in forming landunit (and higher level) averages.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds ! processor bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c, l
    !-----------------------------------------------------------------------

    if (allocated(main_scale_c2l)) deallocate(main_scale_c2l)
    allocate(main_scale_c2l(bounds%begc:bounds%endc, num_c2l_types))

    do c = bounds%begc, bounds%endc
       main_scale_c2l(c,:) = 1.0_r8
       l = col_pp%landunit(c)
       if (lun_pp%urbpoi(l)) then
          if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall) then
             main_scale_c2l(c,urbanf) = 3.0 * lun_pp%canyon_hwr(l)
             main_scale_c2l(c,urbans) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
          else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
             main_scale_c2l(c,urbanf) = 3.0_r8
             main_scale_c2l(c,urbans) = 3.0 / (2.*lun_pp%canyon_hwr(l) + 1.)
          end if
       end if
    end do

    !$acc update device(main_scale_c2l(:,:))

  end subroutine init_scale_c2l

  !-----------------------------------------------------------------------
  integer function c2l_type_id(name)
    ! Integer c2l scale type for a scale type name
    character(len=*), intent(in) :: name
    c2l_type_id = unity
    select case (trim(name))
    case ('unity')
       c2l_type_id = unity
    case ('urbanf')
       c2l_type_id = urbanf
    case ('urbans')
       c2l_type_id = urbans
    case default
       write(iulog,*) 'subgridAveMod error: c2l scale type ', trim(name), ' not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
  end function c2l_type_id

  !-----------------------------------------------------------------------
  integer function l2g_type_id(name)
    ! Integer l2g (or l2t) scale type for a scale type name
    character(len=*), intent(in) :: name
    l2g_type_id = unity
    select case (trim(name))
    case ('unity')
       l2g_type_id = unity
    case ('natveg')
       l2g_type_id = natveg
    case ('veg')
       l2g_type_id = veg
    case ('ice')
       l2g_type_id = ice
    case ('nonurb')
       l2g_type_id = nonurb
    case ('lake')
       l2g_type_id = lake
    case default
       write(iulog,*) 'subgridAveMod error: l2g scale type ', trim(name), ' not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
  end function l2g_type_id

  !-----------------------------------------------------------------------
  integer function p2c_type_id(name)
    ! Integer p2c scale type for a scale type name (only 'unity')
    character(len=*), intent(in) :: name
    if (trim(name) /= 'unity') then
       write(iulog,*) 'subgridAveMod error: p2c scale type ', trim(name), ' not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    p2c_type_id = unity
  end function p2c_type_id

  !-----------------------------------------------------------------------
  subroutine check_scale_types(subname, p2c_scale_type, c2l_scale_type, l2g_scale_type)
    ! Abort on a scale type outside the lookups
    character(len=*), intent(in) :: subname
    integer, intent(in), optional :: p2c_scale_type, c2l_scale_type, l2g_scale_type

    if (present(p2c_scale_type)) then
       if (p2c_scale_type /= unity) call bad_type('p2c', p2c_scale_type)
    end if
    if (present(c2l_scale_type)) then
       if (c2l_scale_type < 1 .or. c2l_scale_type > num_c2l_types) call bad_type('c2l', c2l_scale_type)
    end if
    if (present(l2g_scale_type)) then
       if (l2g_scale_type < 1 .or. l2g_scale_type > num_l2g_types) call bad_type('l2g', l2g_scale_type)
    end if
  contains
    subroutine bad_type(kind, scale_type)
      character(len=*), intent(in) :: kind
      integer, intent(in) :: scale_type
      write(iulog,*) subname, ' error: ', kind, ' scale type ', scale_type, ' not supported'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    end subroutine bad_type
  end subroutine check_scale_types

  !-----------------------------------------------------------------------
  subroutine check_sumwt(subname, err_index, elmlevel)
    ! Abort if an average found a sum of weights greater than 1
    character(len=*), intent(in) :: subname
    integer, intent(in) :: err_index           ! parent with sumwt > 1 (largest), or 0
    character(len=*), intent(in) :: elmlevel   ! level of the parent
    if (err_index /= 0) then
       write(iulog,*) subname, ' error: sumwt is greater than 1.0 at ', trim(elmlevel), ' = ', err_index
       call endrun(decomp_index=err_index, elmlevel=elmlevel, msg=errMsg(__FILE__, __LINE__))
    end if
  end subroutine check_sumwt

  !-----------------------------------------------------------------------
  pure real(r8) function wavg(sum_arr, sumwt)
    !$acc routine seq
    ! Weighted average, spval where nothing contributed
    real(r8), intent(in) :: sum_arr, sumwt
    if (sumwt /= 0._r8) then
       wavg = sum_arr/sumwt
    else
       wavg = spval
    end if
  end function wavg

  !-----------------------------------------------------------------------
  subroutine p2c_1d (bounds, parr, carr, p2c_scale_type)
     ! !DESCRIPTION:
     ! Perfrom subgrid-average from pfts to columns.
     ! Unlike p2c_2d, every patch with nonzero weight contributes (no active or spval check).
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     real(r8), intent(in)  :: parr( bounds%begp: )   ! patch array
     real(r8), intent(out) :: carr( bounds%begc: )   ! column array
     integer , intent(in)  :: p2c_scale_type         ! scale type
     !
     ! !LOCAL VARIABLES:
     integer  :: p,c          ! indices
     integer  :: err_index    ! column with sumwt > 1, or 0
     real(r8) :: sumwt, sum_arr
     !------------------------------------------------------------------------

     SHR_ASSERT_ALL((ubound(parr) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
     SHR_ASSERT_ALL((ubound(carr) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
     call check_scale_types('p2c_1d', p2c_scale_type=p2c_scale_type)

     err_index = 0
     !$acc parallel loop gang vector default(present) &
     !$acc private(sumwt,sum_arr) reduction(max:err_index)
     do c = bounds%begc, bounds%endc
        sumwt = 0._r8
        sum_arr = 0._r8
        !$acc loop seq
        do p = col_pp%pfti(c), col_pp%pftf(c)
           if (veg_pp%wtcol(p) /= 0._r8) then
              sum_arr = sum_arr + parr(p) * veg_pp%wtcol(p)
              sumwt = sumwt + veg_pp%wtcol(p)
           end if
        end do
        if (sumwt > sumwt_tol) err_index = max(err_index, c)
        carr(c) = wavg(sum_arr, sumwt)
     end do
     call check_sumwt('p2c_1d', err_index, namec)

  end subroutine p2c_1d

  !-----------------------------------------------------------------------
  subroutine p2c_2d (bounds, num2d, parr, carr, p2c_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to columns.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds
    integer           , intent(in)  :: num2d                     ! size of second dimension
    real(r8)          , intent(in)  :: parr( bounds%begp: , 1: ) ! patch array
    real(r8)          , intent(out) :: carr( bounds%begc: , 1: ) ! column array
    integer           , intent(in)  :: p2c_scale_type            ! scale type
    !
    ! !LOCAL VARIABLES:
    integer  :: j,p,c        ! indices
    integer  :: err_index    ! column with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(parr) == (/bounds%endp, num2d/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(carr) == (/bounds%endc, num2d/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('p2c_2d', p2c_scale_type=p2c_scale_type)

    err_index = 0
    !$acc parallel loop gang vector collapse(2) default(present) &
    !$acc private(sumwt,sum_arr) reduction(max:err_index)
    do j = 1, num2d
       do c = bounds%begc, bounds%endc
          sumwt = 0._r8
          sum_arr = 0._r8
          !$acc loop seq
          do p = col_pp%pfti(c), col_pp%pftf(c)
             if (veg_pp%active(p) .and. veg_pp%wtcol(p) /= 0._r8) then
                if (parr(p,j) /= spval) then
                   sum_arr = sum_arr + parr(p,j) * veg_pp%wtcol(p)
                   sumwt = sumwt + veg_pp%wtcol(p)
                end if
             end if
          end do
          if (sumwt > sumwt_tol) err_index = max(err_index, c)
          carr(c,j) = wavg(sum_arr, sumwt)
       end do
    end do
    call check_sumwt('p2c_2d', err_index, namec)

  end subroutine p2c_2d

  !-----------------------------------------------------------------------
  subroutine p2c_1d_filter(bounds,numfc, filterc, pftarr, colarr)
    !!
    !! !DESCRIPTION:
    !! perform pft to column averaging for single level pft arrays
    !! (weighted sum over patches with positive weight, no normalization)
    !!
    !! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
    real(r8), intent(in)  :: pftarr(bounds%begp:)
    real(r8), intent(inout) :: colarr(bounds%begc:)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,p  ! indices
    real(r8) :: sum_patch
    !-----------------------------------------------------------------------

    !$acc parallel loop gang vector default(present) &
    !$acc private(c,sum_patch)
    do fc = 1,numfc
       c = filterc(fc)
       sum_patch = 0._r8
       !$acc loop seq
       do p = col_pp%pfti(c), col_pp%pftf(c)
          if (veg_pp%wtcol(p) > 0._r8) sum_patch = sum_patch + pftarr(p) * veg_pp%wtcol(p)
       end do
       colarr(c) = sum_patch
    end do

  end subroutine p2c_1d_filter

  !-----------------------------------------------------------------------
  subroutine p2l_1d (bounds, parr, larr, p2c_scale_type, c2l_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to landunits
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )  ! input pft array
    real(r8), intent(out) :: larr( bounds%begl: )  ! output landunit array
    integer , intent(in)  :: p2c_scale_type        ! scale factor type for averaging
    integer , intent(in)  :: c2l_scale_type        ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,l        ! indices
    integer  :: err_index    ! landunit with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr, scale_c2l
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(parr) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(larr) == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('p2l_1d', p2c_scale_type, c2l_scale_type)

    err_index = 0
    !$acc parallel loop gang vector default(present) &
    !$acc private(sumwt,sum_arr,scale_c2l,c) reduction(max:err_index)
    do l = bounds%begl, bounds%endl
       sumwt = 0._r8
       sum_arr = 0._r8
       !$acc loop seq
       do p = lun_pp%pfti(l), lun_pp%pftf(l)
          if (veg_pp%active(p) .and. veg_pp%wtlunit(p) /= 0._r8) then
             c = veg_pp%column(p)
             scale_c2l = main_scale_c2l(c,c2l_scale_type)
             if (parr(p) /= spval .and. scale_c2l /= spval) then
                sum_arr = sum_arr + parr(p) * scale_c2l * veg_pp%wtlunit(p)
                sumwt = sumwt + veg_pp%wtlunit(p)
             end if
          end if
       end do
       if (sumwt > sumwt_tol) err_index = max(err_index, l)
       larr(l) = wavg(sum_arr, sumwt)
    end do
    call check_sumwt('p2l_1d', err_index, namel)

  end subroutine p2l_1d

  !-----------------------------------------------------------------------
  subroutine p2l_2d(bounds, num2d, parr, larr, p2c_scale_type, c2l_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to landunits
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                      ! size of second dimension
    real(r8), intent(in)  :: parr( bounds%begp: , 1: )  ! input pft array
    real(r8), intent(out) :: larr( bounds%begl: , 1: )  ! output landunit array
    integer , intent(in)  :: p2c_scale_type             ! scale factor type for averaging
    integer , intent(in)  :: c2l_scale_type             ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,p,c,l      ! indices
    integer  :: err_index    ! landunit with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr, scale_c2l
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(parr) == (/bounds%endp, num2d/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(larr) == (/bounds%endl, num2d/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('p2l_2d', p2c_scale_type, c2l_scale_type)

    err_index = 0
    !$acc parallel loop gang vector collapse(2) default(present) &
    !$acc private(sumwt,sum_arr,scale_c2l,c) reduction(max:err_index)
    do j = 1, num2d
       do l = bounds%begl, bounds%endl
          sumwt = 0._r8
          sum_arr = 0._r8
          !$acc loop seq
          do p = lun_pp%pfti(l), lun_pp%pftf(l)
             if (veg_pp%active(p) .and. veg_pp%wtlunit(p) /= 0._r8) then
                c = veg_pp%column(p)
                scale_c2l = main_scale_c2l(c,c2l_scale_type)
                if (parr(p,j) /= spval .and. scale_c2l /= spval) then
                   sum_arr = sum_arr + parr(p,j) * scale_c2l * veg_pp%wtlunit(p)
                   sumwt = sumwt + veg_pp%wtlunit(p)
                end if
             end if
          end do
          if (sumwt > sumwt_tol) err_index = max(err_index, l)
          larr(l,j) = wavg(sum_arr, sumwt)
       end do
    end do
    call check_sumwt('p2l_2d', err_index, namel)

  end subroutine p2l_2d

  !-----------------------------------------------------------------------
  subroutine p2g_1d(bounds, parr, garr, p2c_scale_type, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in), contiguous :: parr( bounds%begp: )  ! input pft array
    real(r8), intent(out), contiguous :: garr( bounds%begg: )  ! output gridcell array
    integer , intent(in)  :: p2c_scale_type        ! scale factor type for averaging
    integer , intent(in)  :: c2l_scale_type        ! scale factor type for averaging
    integer , intent(in)  :: l2g_scale_type        ! scale factor type for averaging
    !
    !  !LOCAL VARIABLES:
    integer  :: chunk, nchunk, g0, g1 ! chunks of gridcells
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(parr) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(garr) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('p2g_1d', p2c_scale_type, c2l_scale_type, l2g_scale_type)

    nchunk = (bounds%endg - bounds%begg) / gridcell_chunk + 1
    !$acc parallel loop gang default(present) private(g0, g1)
    do chunk = 1, nchunk
       g0 = bounds%begg + (chunk-1) * gridcell_chunk
       g1 = min(g0 + gridcell_chunk - 1, bounds%endg)
       if (c2l_scale_type == unity .and. l2g_scale_type == unity) then
          call p2g_range_unity(elm_subgrid_view, g0, g1, parr, bounds%begp, garr(g0:g1))
       else
          call p2g_range(elm_subgrid_view, g0, g1, parr, bounds%begp, c2l_scale_type, l2g_scale_type, garr(g0:g1))
       end if
    end do

  end subroutine p2g_1d

  !-----------------------------------------------------------------------
  subroutine p2g_2d(bounds, num2d, parr, garr, p2c_scale_type, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in), contiguous :: parr( bounds%begp: , 1: ) ! input pft array
    real(r8), intent(out), contiguous :: garr( bounds%begg: , 1: ) ! output gridcell array
    integer , intent(in)  :: p2c_scale_type            ! scale factor type for averaging
    integer , intent(in)  :: c2l_scale_type            ! scale factor type for averaging
    integer , intent(in)  :: l2g_scale_type            ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j            ! level
    integer  :: chunk, nchunk, g0, g1 ! chunks of gridcells
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(parr) == (/bounds%endp, num2d/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(garr) == (/bounds%endg, num2d/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('p2g_2d', p2c_scale_type, c2l_scale_type, l2g_scale_type)

    nchunk = (bounds%endg - bounds%begg) / gridcell_chunk + 1
    !$acc parallel loop gang collapse(2) default(present) private(g0, g1)
    do j = 1, num2d
       do chunk = 1, nchunk
          g0 = bounds%begg + (chunk-1) * gridcell_chunk
          g1 = min(g0 + gridcell_chunk - 1, bounds%endg)
          if (c2l_scale_type == unity .and. l2g_scale_type == unity) then
             call p2g_range_unity(elm_subgrid_view, g0, g1, parr(:,j), bounds%begp, garr(g0:g1,j))
          else
             call p2g_range(elm_subgrid_view, g0, g1, parr(:,j), bounds%begp, c2l_scale_type, l2g_scale_type, garr(g0:g1,j))
          end if
       end do
    end do

  end subroutine p2g_2d

  !-----------------------------------------------------------------------
  subroutine c2l_1d (bounds, carr, larr, c2l_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from columns to landunits
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: carr( bounds%begc: )  ! input column array
    real(r8), intent(out) :: larr( bounds%begl: )  ! output landunit array
    integer , intent(in)  :: c2l_scale_type        ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l          ! indices
    integer  :: err_index    ! landunit with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr, scale_c2l
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(carr) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(larr) == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('c2l_1d', c2l_scale_type=c2l_scale_type)

    err_index = 0
    !$acc parallel loop gang vector default(present) &
    !$acc private(sumwt,sum_arr,scale_c2l) reduction(max:err_index)
    do l = bounds%begl, bounds%endl
       sumwt = 0._r8
       sum_arr = 0._r8
       !$acc loop seq
       do c = lun_pp%coli(l), lun_pp%colf(l)
          if (col_pp%active(c) .and. col_pp%wtlunit(c) /= 0._r8) then
             scale_c2l = main_scale_c2l(c,c2l_scale_type)
             if (carr(c) /= spval .and. scale_c2l /= spval) then
                sum_arr = sum_arr + carr(c) * scale_c2l * col_pp%wtlunit(c)
                sumwt = sumwt + col_pp%wtlunit(c)
             end if
          end if
       end do
       if (sumwt > sumwt_tol) err_index = max(err_index, l)
       larr(l) = wavg(sum_arr, sumwt)
    end do
    call check_sumwt('c2l_1d', err_index, namel)

  end subroutine c2l_1d

  !-----------------------------------------------------------------------
  subroutine c2l_2d (bounds, num2d, carr, larr, c2l_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from columns to landunits
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: carr( bounds%begc: , 1: ) ! input column array
    real(r8), intent(out) :: larr( bounds%begl: , 1: ) ! output landunit array
    integer , intent(in)  :: c2l_scale_type            ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c        ! indices
    integer  :: err_index    ! landunit with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr, scale_c2l
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(carr) == (/bounds%endc, num2d/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(larr) == (/bounds%endl, num2d/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('c2l_2d', c2l_scale_type=c2l_scale_type)

    err_index = 0
    !$acc parallel loop gang vector collapse(2) default(present) &
    !$acc private(sumwt,sum_arr,scale_c2l) reduction(max:err_index)
    do j = 1, num2d
       do l = bounds%begl, bounds%endl
          sumwt = 0._r8
          sum_arr = 0._r8
          !$acc loop seq
          do c = lun_pp%coli(l), lun_pp%colf(l)
             if (col_pp%active(c) .and. col_pp%wtlunit(c) /= 0._r8) then
                scale_c2l = main_scale_c2l(c,c2l_scale_type)
                if (carr(c,j) /= spval .and. scale_c2l /= spval) then
                   sum_arr = sum_arr + carr(c,j) * scale_c2l * col_pp%wtlunit(c)
                   sumwt = sumwt + col_pp%wtlunit(c)
                end if
             end if
          end do
          if (sumwt > sumwt_tol) err_index = max(err_index, l)
          larr(l,j) = wavg(sum_arr, sumwt)
       end do
    end do
    call check_sumwt('c2l_2d', err_index, namel)

  end subroutine c2l_2d

  !-----------------------------------------------------------------------
  subroutine c2g_1d(bounds, carr, garr, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from columns to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in), contiguous :: carr( bounds%begc: )  ! input column array
    real(r8), intent(out), contiguous :: garr( bounds%begg: )  ! output gridcell array
    integer , intent(in)  :: c2l_scale_type        ! scale factor type for averaging
    integer , intent(in)  :: l2g_scale_type        ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: chunk, nchunk, g0, g1 ! chunks of gridcells
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(carr) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(garr) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('c2g_1d', c2l_scale_type=c2l_scale_type, l2g_scale_type=l2g_scale_type)

    nchunk = (bounds%endg - bounds%begg) / gridcell_chunk + 1
    !$acc parallel loop gang default(present) private(g0, g1)
    do chunk = 1, nchunk
       g0 = bounds%begg + (chunk-1) * gridcell_chunk
       g1 = min(g0 + gridcell_chunk - 1, bounds%endg)
       if (c2l_scale_type == unity .and. l2g_scale_type == unity) then
          call c2g_range_unity(elm_subgrid_view, g0, g1, carr, bounds%begc, garr(g0:g1))
       else
          call c2g_range(elm_subgrid_view, g0, g1, carr, bounds%begc, c2l_scale_type, l2g_scale_type, garr(g0:g1))
       end if
    end do

  end subroutine c2g_1d

  !-----------------------------------------------------------------------
  subroutine c2g_2d(bounds, num2d, carr, garr, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from columns to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in), contiguous :: carr( bounds%begc: , 1: ) ! input column array
    real(r8), intent(out), contiguous :: garr( bounds%begg: , 1: ) ! output gridcell array
    integer , intent(in)  :: c2l_scale_type            ! scale factor type for averaging
    integer , intent(in)  :: l2g_scale_type            ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j            ! level
    integer  :: chunk, nchunk, g0, g1 ! chunks of gridcells
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(carr) == (/bounds%endc, num2d/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(garr) == (/bounds%endg, num2d/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('c2g_2d', c2l_scale_type=c2l_scale_type, l2g_scale_type=l2g_scale_type)

    nchunk = (bounds%endg - bounds%begg) / gridcell_chunk + 1
    !$acc parallel loop gang collapse(2) default(present) private(g0, g1)
    do j = 1, num2d
       do chunk = 1, nchunk
          g0 = bounds%begg + (chunk-1) * gridcell_chunk
          g1 = min(g0 + gridcell_chunk - 1, bounds%endg)
          if (c2l_scale_type == unity .and. l2g_scale_type == unity) then
             call c2g_range_unity(elm_subgrid_view, g0, g1, carr(:,j), bounds%begc, garr(g0:g1,j))
          else
             call c2g_range(elm_subgrid_view, g0, g1, carr(:,j), bounds%begc, c2l_scale_type, l2g_scale_type, garr(g0:g1,j))
          end if
       end do
    end do

  end subroutine c2g_2d

  !-----------------------------------------------------------------------
  subroutine l2g_1d(bounds, larr, garr, l2g_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from landunits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in), contiguous :: larr( bounds%begl: )  ! input landunit array
    real(r8), intent(out), contiguous :: garr( bounds%begg: )  ! output gridcell array
    integer , intent(in)  :: l2g_scale_type        ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: chunk, nchunk, g0, g1 ! chunks of gridcells
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(larr) == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(garr) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('l2g_1d', l2g_scale_type=l2g_scale_type)

    nchunk = (bounds%endg - bounds%begg) / gridcell_chunk + 1
    !$acc parallel loop gang default(present) private(g0, g1)
    do chunk = 1, nchunk
       g0 = bounds%begg + (chunk-1) * gridcell_chunk
       g1 = min(g0 + gridcell_chunk - 1, bounds%endg)
       call l2g_range(elm_subgrid_view, g0, g1, larr, bounds%begl, l2g_scale_type, garr(g0:g1))
    end do

  end subroutine l2g_1d

  !-----------------------------------------------------------------------
  subroutine l2g_2d(bounds, num2d, larr, garr, l2g_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from landunits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in), contiguous :: larr( bounds%begl: , 1: ) ! input landunit array
    real(r8), intent(out), contiguous :: garr( bounds%begg: , 1: ) ! output gridcell array
    integer , intent(in)  :: l2g_scale_type            ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j            ! level
    integer  :: chunk, nchunk, g0, g1 ! chunks of gridcells
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(larr) == (/bounds%endl, num2d/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(garr) == (/bounds%endg, num2d/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('l2g_2d', l2g_scale_type=l2g_scale_type)

    nchunk = (bounds%endg - bounds%begg) / gridcell_chunk + 1
    !$acc parallel loop gang collapse(2) default(present) private(g0, g1)
    do j = 1, num2d
       do chunk = 1, nchunk
          g0 = bounds%begg + (chunk-1) * gridcell_chunk
          g1 = min(g0 + gridcell_chunk - 1, bounds%endg)
          call l2g_range(elm_subgrid_view, g0, g1, larr(:,j), bounds%begl, l2g_scale_type, garr(g0:g1,j))
       end do
    end do

  end subroutine l2g_2d

  !-----------------------------------------------------------------------
  subroutine p2t_1d(bounds, parr, tarr, p2c_scale_type, c2l_scale_type, l2t_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to topounits.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )  ! input pft array
    real(r8), intent(out) :: tarr( bounds%begt: )  ! output topounit array
    integer , intent(in)  :: p2c_scale_type        ! scale factor type for averaging
    integer , intent(in)  :: c2l_scale_type        ! scale factor type for averaging
    integer , intent(in)  :: l2t_scale_type        ! scale factor type for averaging
    !
    !  !LOCAL VARIABLES:
    integer  :: k,p,c,l,t    ! indices
    integer  :: err_index    ! topounit with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr, scale_c2l, scale_l2t
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(parr) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tarr) == (/bounds%endt/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('p2t_1d', p2c_scale_type, c2l_scale_type, l2t_scale_type)

    err_index = 0
    !$acc parallel loop gang vector default(present) &
    !$acc private(sumwt,sum_arr,scale_c2l,scale_l2t,p,c,l) reduction(max:err_index)
    do t = bounds%begt, bounds%endt
       sumwt = 0._r8
       sum_arr = 0._r8
       !$acc loop seq
       do k = top_pft_ptr(t), top_pft_ptr(t+1)-1
          p = top_pft_idx(k)
          if (veg_pp%active(p) .and. veg_pp%wttopounit(p) /= 0._r8) then
             c = veg_pp%column(p)
             l = veg_pp%landunit(p)
             scale_c2l = main_scale_c2l(c,c2l_scale_type)
             scale_l2t = main_scale_l2g_lookup(l2t_scale_type,lun_pp%itype(l))
             if (parr(p) /= spval .and. scale_c2l /= spval .and. scale_l2t /= spval) then
                sum_arr = sum_arr + parr(p) * scale_c2l * scale_l2t * veg_pp%wttopounit(p)
                sumwt = sumwt + veg_pp%wttopounit(p)
             end if
          end if
       end do
       if (sumwt > sumwt_tol) err_index = max(err_index, t)
       tarr(t) = wavg(sum_arr, sumwt)
    end do
    call check_sumwt('p2t_1d', err_index, namet)

  end subroutine p2t_1d

  !-----------------------------------------------------------------------
  subroutine p2t_2d(bounds, num2d, parr, tarr, p2c_scale_type, c2l_scale_type, l2t_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to topounits.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: parr( bounds%begp: , 1: ) ! input pft array
    real(r8), intent(out) :: tarr( bounds%begt: , 1: ) ! output topounit array
    integer , intent(in)  :: p2c_scale_type            ! scale factor type for averaging
    integer , intent(in)  :: c2l_scale_type            ! scale factor type for averaging
    integer , intent(in)  :: l2t_scale_type            ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,k,p,c,l,t  ! indices
    integer  :: err_index    ! topounit with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr, scale_c2l, scale_l2t
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(parr) == (/bounds%endp, num2d/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tarr) == (/bounds%endt, num2d/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('p2t_2d', p2c_scale_type, c2l_scale_type, l2t_scale_type)

    err_index = 0
    !$acc parallel loop gang vector collapse(2) default(present) &
    !$acc private(sumwt,sum_arr,scale_c2l,scale_l2t,p,c,l) reduction(max:err_index)
    do j = 1, num2d
       do t = bounds%begt, bounds%endt
          sumwt = 0._r8
          sum_arr = 0._r8
          !$acc loop seq
          do k = top_pft_ptr(t), top_pft_ptr(t+1)-1
             p = top_pft_idx(k)
             if (veg_pp%active(p) .and. veg_pp%wttopounit(p) /= 0._r8) then
                c = veg_pp%column(p)
                l = veg_pp%landunit(p)
                scale_c2l = main_scale_c2l(c,c2l_scale_type)
                scale_l2t = main_scale_l2g_lookup(l2t_scale_type,lun_pp%itype(l))
                if (parr(p,j) /= spval .and. scale_c2l /= spval .and. scale_l2t /= spval) then
                   sum_arr = sum_arr + parr(p,j) * scale_c2l * scale_l2t * veg_pp%wttopounit(p)
                   sumwt = sumwt + veg_pp%wttopounit(p)
                end if
             end if
          end do
          if (sumwt > sumwt_tol) err_index = max(err_index, t)
          tarr(t,j) = wavg(sum_arr, sumwt)
       end do
    end do
    call check_sumwt('p2t_2d', err_index, namet)

  end subroutine p2t_2d

  !-----------------------------------------------------------------------
  subroutine c2t_1d(bounds, carr, tarr, c2l_scale_type, l2t_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from columns to topounits.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: carr( bounds%begc: )  ! input column array
    real(r8), intent(out) :: tarr( bounds%begt: )  ! output topounit array
    integer , intent(in)  :: c2l_scale_type        ! scale factor type for averaging
    integer , intent(in)  :: l2t_scale_type        ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: k,c,l,t      ! indices
    integer  :: err_index    ! topounit with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr, scale_c2l, scale_l2t
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(carr) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tarr) == (/bounds%endt/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('c2t_1d', c2l_scale_type=c2l_scale_type, l2g_scale_type=l2t_scale_type)

    err_index = 0
    !$acc parallel loop gang vector default(present) &
    !$acc private(sumwt,sum_arr,scale_c2l,scale_l2t,c,l) reduction(max:err_index)
    do t = bounds%begt, bounds%endt
       sumwt = 0._r8
       sum_arr = 0._r8
       !$acc loop seq
       do k = top_col_ptr(t), top_col_ptr(t+1)-1
          c = top_col_idx(k)
          if (col_pp%active(c) .and. col_pp%wttopounit(c) /= 0._r8) then
             l = col_pp%landunit(c)
             scale_c2l = main_scale_c2l(c,c2l_scale_type)
             scale_l2t = main_scale_l2g_lookup(l2t_scale_type,lun_pp%itype(l))
             if (carr(c) /= spval .and. scale_c2l /= spval .and. scale_l2t /= spval) then
                sum_arr = sum_arr + carr(c) * scale_c2l * scale_l2t * col_pp%wttopounit(c)
                sumwt = sumwt + col_pp%wttopounit(c)
             end if
          end if
       end do
       if (sumwt > sumwt_tol) err_index = max(err_index, t)
       tarr(t) = wavg(sum_arr, sumwt)
    end do
    call check_sumwt('c2t_1d', err_index, namet)

  end subroutine c2t_1d

  !-----------------------------------------------------------------------
  subroutine c2t_2d(bounds, num2d, carr, tarr, c2l_scale_type, l2t_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from columns to topounits.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: carr( bounds%begc: , 1: ) ! input column array
    real(r8), intent(out) :: tarr( bounds%begt: , 1: ) ! output topounit array
    integer , intent(in)  :: c2l_scale_type            ! scale factor type for averaging
    integer , intent(in)  :: l2t_scale_type            ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,k,c,l,t    ! indices
    integer  :: err_index    ! topounit with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr, scale_c2l, scale_l2t
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(carr) == (/bounds%endc, num2d/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tarr) == (/bounds%endt, num2d/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('c2t_2d', c2l_scale_type=c2l_scale_type, l2g_scale_type=l2t_scale_type)

    err_index = 0
    !$acc parallel loop gang vector collapse(2) default(present) &
    !$acc private(sumwt,sum_arr,scale_c2l,scale_l2t,c,l) reduction(max:err_index)
    do j = 1, num2d
       do t = bounds%begt, bounds%endt
          sumwt = 0._r8
          sum_arr = 0._r8
          !$acc loop seq
          do k = top_col_ptr(t), top_col_ptr(t+1)-1
             c = top_col_idx(k)
             if (col_pp%active(c) .and. col_pp%wttopounit(c) /= 0._r8) then
                l = col_pp%landunit(c)
                scale_c2l = main_scale_c2l(c,c2l_scale_type)
                scale_l2t = main_scale_l2g_lookup(l2t_scale_type,lun_pp%itype(l))
                if (carr(c,j) /= spval .and. scale_c2l /= spval .and. scale_l2t /= spval) then
                   sum_arr = sum_arr + carr(c,j) * scale_c2l * scale_l2t * col_pp%wttopounit(c)
                   sumwt = sumwt + col_pp%wttopounit(c)
                end if
             end if
          end do
          if (sumwt > sumwt_tol) err_index = max(err_index, t)
          tarr(t,j) = wavg(sum_arr, sumwt)
       end do
    end do
    call check_sumwt('c2t_2d', err_index, namet)

  end subroutine c2t_2d

  !-----------------------------------------------------------------------
  subroutine l2t_1d(bounds, larr, tarr, l2t_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from landunits to topounits.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: larr( bounds%begl: )  ! input landunit array
    real(r8), intent(out) :: tarr( bounds%begt: )  ! output topounit array
    integer , intent(in)  :: l2t_scale_type        ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: k,l,t        ! indices
    integer  :: err_index    ! topounit with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr, scale_l2t
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(larr) == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tarr) == (/bounds%endt/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('l2t_1d', l2g_scale_type=l2t_scale_type)

    err_index = 0
    !$acc parallel loop gang vector default(present) &
    !$acc private(sumwt,sum_arr,scale_l2t,l) reduction(max:err_index)
    do t = bounds%begt, bounds%endt
       sumwt = 0._r8
       sum_arr = 0._r8
       !$acc loop seq
       do k = top_lun_ptr(t), top_lun_ptr(t+1)-1
          l = top_lun_idx(k)
          if (lun_pp%active(l) .and. lun_pp%wttopounit(l) /= 0._r8) then
             scale_l2t = main_scale_l2g_lookup(l2t_scale_type,lun_pp%itype(l))
             if (larr(l) /= spval .and. scale_l2t /= spval) then
                sum_arr = sum_arr + larr(l) * scale_l2t * lun_pp%wttopounit(l)
                sumwt = sumwt + lun_pp%wttopounit(l)
             end if
          end if
       end do
       if (sumwt > sumwt_tol) err_index = max(err_index, t)
       tarr(t) = wavg(sum_arr, sumwt)
    end do
    call check_sumwt('l2t_1d', err_index, namet)

  end subroutine l2t_1d

  !-----------------------------------------------------------------------
  subroutine l2t_2d(bounds, num2d, larr, tarr, l2t_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from landunits to topounits.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: larr( bounds%begl: , 1: ) ! input landunit array
    real(r8), intent(out) :: tarr( bounds%begt: , 1: ) ! output topounit array
    integer , intent(in)  :: l2t_scale_type            ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,k,l,t      ! indices
    integer  :: err_index    ! topounit with sumwt > 1, or 0
    real(r8) :: sumwt, sum_arr, scale_l2t
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(larr) == (/bounds%endl, num2d/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tarr) == (/bounds%endt, num2d/)), errMsg(__FILE__, __LINE__))
    call check_scale_types('l2t_2d', l2g_scale_type=l2t_scale_type)

    err_index = 0
    !$acc parallel loop gang vector collapse(2) default(present) &
    !$acc private(sumwt,sum_arr,scale_l2t,l) reduction(max:err_index)
    do j = 1, num2d
       do t = bounds%begt, bounds%endt
          sumwt = 0._r8
          sum_arr = 0._r8
          !$acc loop seq
          do k = top_lun_ptr(t), top_lun_ptr(t+1)-1
             l = top_lun_idx(k)
             if (lun_pp%active(l) .and. lun_pp%wttopounit(l) /= 0._r8) then
                scale_l2t = main_scale_l2g_lookup(l2t_scale_type,lun_pp%itype(l))
                if (larr(l,j) /= spval .and. scale_l2t /= spval) then
                   sum_arr = sum_arr + larr(l,j) * scale_l2t * lun_pp%wttopounit(l)
                   sumwt = sumwt + lun_pp%wttopounit(l)
                end if
             end if
          end do
          if (sumwt > sumwt_tol) err_index = max(err_index, t)
          tarr(t,j) = wavg(sum_arr, sumwt)
       end do
    end do
    call check_sumwt('l2t_2d', err_index, namet)

  end subroutine l2t_2d

  !-----------------------------------------------------------------------
  subroutine t2g_1d(bounds, tarr, garr, t2g_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from topounits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in), contiguous :: tarr( bounds%begt: )  ! input topounit array
    real(r8), intent(out), contiguous :: garr( bounds%begg: )  ! output gridcell array
    character(len=*), intent(in) :: t2g_scale_type ! scale factor type for averaging (only 'unity')
    !
    ! !LOCAL VARIABLES:
    integer  :: chunk, nchunk, g0, g1 ! chunks of gridcells
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(tarr) == (/bounds%endt/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(garr) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))
    if (trim(t2g_scale_type) /= 'unity') then
       write(iulog,*) 't2g_1d error: scale type ', trim(t2g_scale_type), ' not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    nchunk = (bounds%endg - bounds%begg) / gridcell_chunk + 1
    !$acc parallel loop gang default(present) private(g0, g1)
    do chunk = 1, nchunk
       g0 = bounds%begg + (chunk-1) * gridcell_chunk
       g1 = min(g0 + gridcell_chunk - 1, bounds%endg)
       call t2g_range(elm_subgrid_view, g0, g1, tarr, bounds%begt, garr(g0:g1))
    end do

  end subroutine t2g_1d

  !-----------------------------------------------------------------------
  subroutine t2g_2d(bounds, num2d, tarr, garr, t2g_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from topounits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in), contiguous :: tarr( bounds%begt: , 1: ) ! input topounit array
    real(r8), intent(out), contiguous :: garr( bounds%begg: , 1: ) ! output gridcell array
    character(len=*), intent(in) :: t2g_scale_type     ! scale factor type for averaging (only 'unity')
    !
    ! !LOCAL VARIABLES:
    integer  :: j            ! level
    integer  :: chunk, nchunk, g0, g1 ! chunks of gridcells
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(tarr) == (/bounds%endt, num2d/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(garr) == (/bounds%endg, num2d/)), errMsg(__FILE__, __LINE__))
    if (trim(t2g_scale_type) /= 'unity') then
       write(iulog,*) 't2g_2d error: scale type ', trim(t2g_scale_type), ' not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    nchunk = (bounds%endg - bounds%begg) / gridcell_chunk + 1
    !$acc parallel loop gang collapse(2) default(present) private(g0, g1)
    do j = 1, num2d
       do chunk = 1, nchunk
          g0 = bounds%begg + (chunk-1) * gridcell_chunk
          g1 = min(g0 + gridcell_chunk - 1, bounds%endg)
          call t2g_range(elm_subgrid_view, g0, g1, tarr(:,j), bounds%begt, garr(g0:g1,j))
       end do
    end do

  end subroutine t2g_2d

  !-----------------------------------------------------------------------
  ! Versions taking the scale type names
  !-----------------------------------------------------------------------
  subroutine p2c_1d_name(bounds, parr, carr, p2c_scale_type)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )
    real(r8), intent(out) :: carr( bounds%begc: )
    character(len=*), intent(in) :: p2c_scale_type
    call p2c_1d(bounds, parr, carr, p2c_type_id(p2c_scale_type))
  end subroutine p2c_1d_name

  subroutine p2c_2d_name(bounds, num2d, parr, carr, p2c_scale_type)
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d
    real(r8), intent(in)  :: parr( bounds%begp: , 1: )
    real(r8), intent(out) :: carr( bounds%begc: , 1: )
    character(len=*), intent(in) :: p2c_scale_type
    call p2c_2d(bounds, num2d, parr, carr, p2c_type_id(p2c_scale_type))
  end subroutine p2c_2d_name

  subroutine p2l_1d_name(bounds, parr, larr, p2c_scale_type, c2l_scale_type)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )
    real(r8), intent(out) :: larr( bounds%begl: )
    character(len=*), intent(in) :: p2c_scale_type, c2l_scale_type
    call p2l_1d(bounds, parr, larr, p2c_type_id(p2c_scale_type), c2l_type_id(c2l_scale_type))
  end subroutine p2l_1d_name

  subroutine p2l_2d_name(bounds, num2d, parr, larr, p2c_scale_type, c2l_scale_type)
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d
    real(r8), intent(in)  :: parr( bounds%begp: , 1: )
    real(r8), intent(out) :: larr( bounds%begl: , 1: )
    character(len=*), intent(in) :: p2c_scale_type, c2l_scale_type
    call p2l_2d(bounds, num2d, parr, larr, p2c_type_id(p2c_scale_type), c2l_type_id(c2l_scale_type))
  end subroutine p2l_2d_name

  subroutine p2g_1d_name(bounds, parr, garr, p2c_scale_type, c2l_scale_type, l2g_scale_type)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )
    real(r8), intent(out) :: garr( bounds%begg: )
    character(len=*), intent(in) :: p2c_scale_type, c2l_scale_type, l2g_scale_type
    call p2g_1d(bounds, parr, garr, p2c_type_id(p2c_scale_type), c2l_type_id(c2l_scale_type), &
         l2g_type_id(l2g_scale_type))
  end subroutine p2g_1d_name

  subroutine p2g_2d_name(bounds, num2d, parr, garr, p2c_scale_type, c2l_scale_type, l2g_scale_type)
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d
    real(r8), intent(in)  :: parr( bounds%begp: , 1: )
    real(r8), intent(out) :: garr( bounds%begg: , 1: )
    character(len=*), intent(in) :: p2c_scale_type, c2l_scale_type, l2g_scale_type
    call p2g_2d(bounds, num2d, parr, garr, p2c_type_id(p2c_scale_type), c2l_type_id(c2l_scale_type), &
         l2g_type_id(l2g_scale_type))
  end subroutine p2g_2d_name

  subroutine c2l_1d_name(bounds, carr, larr, c2l_scale_type)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: carr( bounds%begc: )
    real(r8), intent(out) :: larr( bounds%begl: )
    character(len=*), intent(in) :: c2l_scale_type
    call c2l_1d(bounds, carr, larr, c2l_type_id(c2l_scale_type))
  end subroutine c2l_1d_name

  subroutine c2l_2d_name(bounds, num2d, carr, larr, c2l_scale_type)
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d
    real(r8), intent(in)  :: carr( bounds%begc: , 1: )
    real(r8), intent(out) :: larr( bounds%begl: , 1: )
    character(len=*), intent(in) :: c2l_scale_type
    call c2l_2d(bounds, num2d, carr, larr, c2l_type_id(c2l_scale_type))
  end subroutine c2l_2d_name

  subroutine c2g_1d_name(bounds, carr, garr, c2l_scale_type, l2g_scale_type)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: carr( bounds%begc: )
    real(r8), intent(out) :: garr( bounds%begg: )
    character(len=*), intent(in) :: c2l_scale_type, l2g_scale_type
    call c2g_1d(bounds, carr, garr, c2l_type_id(c2l_scale_type), l2g_type_id(l2g_scale_type))
  end subroutine c2g_1d_name

  subroutine c2g_2d_name(bounds, num2d, carr, garr, c2l_scale_type, l2g_scale_type)
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d
    real(r8), intent(in)  :: carr( bounds%begc: , 1: )
    real(r8), intent(out) :: garr( bounds%begg: , 1: )
    character(len=*), intent(in) :: c2l_scale_type, l2g_scale_type
    call c2g_2d(bounds, num2d, carr, garr, c2l_type_id(c2l_scale_type), l2g_type_id(l2g_scale_type))
  end subroutine c2g_2d_name

  subroutine l2g_1d_name(bounds, larr, garr, l2g_scale_type)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: larr( bounds%begl: )
    real(r8), intent(out) :: garr( bounds%begg: )
    character(len=*), intent(in) :: l2g_scale_type
    call l2g_1d(bounds, larr, garr, l2g_type_id(l2g_scale_type))
  end subroutine l2g_1d_name

  subroutine l2g_2d_name(bounds, num2d, larr, garr, l2g_scale_type)
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d
    real(r8), intent(in)  :: larr( bounds%begl: , 1: )
    real(r8), intent(out) :: garr( bounds%begg: , 1: )
    character(len=*), intent(in) :: l2g_scale_type
    call l2g_2d(bounds, num2d, larr, garr, l2g_type_id(l2g_scale_type))
  end subroutine l2g_2d_name

  subroutine p2t_1d_name(bounds, parr, tarr, p2c_scale_type, c2l_scale_type, l2t_scale_type)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )
    real(r8), intent(out) :: tarr( bounds%begt: )
    character(len=*), intent(in) :: p2c_scale_type, c2l_scale_type, l2t_scale_type
    call p2t_1d(bounds, parr, tarr, p2c_type_id(p2c_scale_type), c2l_type_id(c2l_scale_type), &
         l2g_type_id(l2t_scale_type))
  end subroutine p2t_1d_name

  subroutine p2t_2d_name(bounds, num2d, parr, tarr, p2c_scale_type, c2l_scale_type, l2t_scale_type)
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d
    real(r8), intent(in)  :: parr( bounds%begp: , 1: )
    real(r8), intent(out) :: tarr( bounds%begt: , 1: )
    character(len=*), intent(in) :: p2c_scale_type, c2l_scale_type, l2t_scale_type
    call p2t_2d(bounds, num2d, parr, tarr, p2c_type_id(p2c_scale_type), c2l_type_id(c2l_scale_type), &
         l2g_type_id(l2t_scale_type))
  end subroutine p2t_2d_name

  subroutine c2t_1d_name(bounds, carr, tarr, c2l_scale_type, l2t_scale_type)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: carr( bounds%begc: )
    real(r8), intent(out) :: tarr( bounds%begt: )
    character(len=*), intent(in) :: c2l_scale_type, l2t_scale_type
    call c2t_1d(bounds, carr, tarr, c2l_type_id(c2l_scale_type), l2g_type_id(l2t_scale_type))
  end subroutine c2t_1d_name

  subroutine c2t_2d_name(bounds, num2d, carr, tarr, c2l_scale_type, l2t_scale_type)
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d
    real(r8), intent(in)  :: carr( bounds%begc: , 1: )
    real(r8), intent(out) :: tarr( bounds%begt: , 1: )
    character(len=*), intent(in) :: c2l_scale_type, l2t_scale_type
    call c2t_2d(bounds, num2d, carr, tarr, c2l_type_id(c2l_scale_type), l2g_type_id(l2t_scale_type))
  end subroutine c2t_2d_name

  subroutine l2t_1d_name(bounds, larr, tarr, l2t_scale_type)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: larr( bounds%begl: )
    real(r8), intent(out) :: tarr( bounds%begt: )
    character(len=*), intent(in) :: l2t_scale_type
    call l2t_1d(bounds, larr, tarr, l2g_type_id(l2t_scale_type))
  end subroutine l2t_1d_name

  subroutine l2t_2d_name(bounds, num2d, larr, tarr, l2t_scale_type)
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d
    real(r8), intent(in)  :: larr( bounds%begl: , 1: )
    real(r8), intent(out) :: tarr( bounds%begt: , 1: )
    character(len=*), intent(in) :: l2t_scale_type
    call l2t_2d(bounds, num2d, larr, tarr, l2g_type_id(l2t_scale_type))
  end subroutine l2t_2d_name

end module subgridAveMod
