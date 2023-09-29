module subgridAveMod

#include "shr_assert.h"
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Utilities to perfrom subgrid averaging
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use column_varcon , only : icol_roof, icol_sunwall, icol_shadewall
  use column_varcon , only : icol_road_perv , icol_road_imperv
  use elm_varcon    , only : grlnd, nameg, namet, namel, namec, namep,spval 
  use elm_varctl    , only : iulog
  use decompMod     , only : bounds_type
  use TopounitType  , only : top_pp
  use LandunitType  , only : lun_pp
  use ColumnType    , only : col_pp,column_physical_properties
  use VegetationType, only : veg_pp
  use abortutils    , only : endrun
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  save

  integer, public, parameter :: unity = 1, natveg = 2, veg = 3
  integer, public, parameter :: ice = 4, nonurb = 5, lake = 6
  integer, public, parameter :: urbanf = 2, urbans = 3
  !
  ! !PUBLIC MEMBER FUNCTIONS:
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

  !! Intefaces for openacc compatible versions of subroutines were created
  !! may consolidate once ELM port is finished
  interface p2c
     module procedure p2c_1d
     module procedure p2c_1d_gpu
     module procedure p2c_2d
   !   module procedure p2c_2d_gpu
     module procedure p2c_1d_filter
     module procedure p2c_2d_filter
  end interface
  interface p2l
     module procedure p2l_1d
     module procedure p2l_2d
  end interface
  interface p2g
     module procedure p2g_1d
     module procedure p2g_1d_gpu
     module procedure p2g_2d
     module procedure p2g_2d_gpu
  end interface
  interface c2l
     module procedure c2l_1d
     module procedure c2l_2d
  end interface
  interface c2g
     module procedure c2g_1d
     module procedure c2g_2d
     module procedure c2g_1d_gpu
     module procedure c2g_2d_gpu
  end interface
  interface l2g
     module procedure l2g_1d
     module procedure l2g_2d
     module procedure l2g_1d_gpu
     module procedure l2g_2d_gpu
  end interface
  interface p2t
     module procedure p2t_1d
     module procedure p2t_2d
  end interface
  interface c2t
     module procedure c2t_1d
     module procedure c2t_2d
  end interface
  interface l2t
     module procedure l2t_1d
     module procedure l2t_2d
  end interface
  interface t2g
     module procedure t2g_1d
     module procedure t2g_2d
     module procedure t2g_1d_gpu
     module procedure t2g_2d_gpu
  end interface


  ! Created this interface to avoid repeated code
  ! implemented similar to create_scale_l2g
  interface create_scale_c2l
    module procedure create_scale_c2l
    module procedure create_scale_c2l_gpu
  end interface

  public :: p2c_1d_filter_parallel
  public :: p2c_1d_parallel 
  public :: p2c_2d_parallel
  public :: c2g_1d_parallel
  !
  public :: initialize_scale_l2g_lookup
  public :: initialize_scale_c2l
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: build_scale_l2g
  private :: create_scale_l2g_lookup
  private :: build_scale_l2t
  private :: create_scale_l2t_lookup

  ! New arrays created to stop unnecessary calls for fixed scaling quantities
  real(r8), allocatable, private :: main_scale_l2g_lookup(:,:) !dimensions (scale_type, lun_pp%itype)
  !$acc declare create(main_scale_l2g_lookup(:,:))
  real(r8), allocatable, private :: main_scale_c2l(:,:) ! dimensions = (col, scale_type) 
  !$acc declare create(main_scale_c2l(:,:))
  
  ! WJS (10-14-11): TODO:
  !
  ! - I believe that scale_p2c, scale_c2l and scale_l2g should be included in the sumwt
  ! accumulations (e.g., sumwt = sumwt + wtgcell * scale_p2c * scale_c2l * scale_l2g), but
  ! that requires some more thought to (1) make sure that is correct, and (2) make sure it
  ! doesn't break the urban scaling. (See also my notes in create_scale_l2g_lookup.)
  !   - Once that is done, you could use a scale of 0, avoiding the need for the use of
  !   spval and the special checks that requires.
  !
  ! - Currently, there is a lot of repeated code to calculate scale_c2l. This should be
  ! cleaned up.
  !   - At a minimum, should collect the repeated code into a subroutine to eliminate this
  !   repitition
  !   - The best thing might be to use a lookup array, as is done for scale_l2g
  ! -----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine p2c_1d (bounds, parr, carr, p2c_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to columns.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
      !$acc routine seq 
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )         ! patch array
    real(r8), intent(out) :: carr( bounds%begc: )         ! column array
    character(len=*), intent(in) :: p2c_scale_type ! scale type
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,index                       ! indices
    logical  :: found                           ! temporary for error check
    real(r8) :: sumwt(bounds%begc:bounds%endc)  ! sum of weights
    !------------------------------------------------------------------------

    carr(bounds%begc:bounds%endc) = spval
    sumwt(bounds%begc:bounds%endc) = 0._r8
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .and. veg_pp%wtcol(p) /= 0._r8) then
          if (parr(p) /= spval) then
             c = veg_pp%column(p)
             if (sumwt(c) == 0._r8) carr(c) = 0._r8
             carr(c) = carr(c) + parr(p) * veg_pp%wtcol(p)
             sumwt(c) = sumwt(c) + veg_pp%wtcol(p)
          end if
       end if
    end do
    
    found = .false.
    do c = bounds%begc,bounds%endc
       if (sumwt(c) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = c
       else if (sumwt(c) /= 0._r8) then
          carr(c) = carr(c)/sumwt(c)
       end if
    end do

  end subroutine p2c_1d

  !-----------------------------------------------------------------------
  subroutine p2c_1d_gpu (bounds, parr, carr, p2c_scale_type)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to columns.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )         ! patch array
    real(r8), intent(out) :: carr( bounds%begc: )         ! column array
    integer, intent(in) :: p2c_scale_type ! scale type
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,index                       ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: sumwt(bounds%begc:bounds%endc)     ! sum of weights
    !------------------------------------------------------------------------

    carr(bounds%begc:bounds%endc) = spval
    sumwt(bounds%begc:bounds%endc) = 0._r8
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .and. veg_pp%wtcol(p) /= 0._r8) then
          if (parr(p) /= spval) then
             c = veg_pp%column(p)
             if (sumwt(c) == 0._r8) carr(c) = 0._r8
             carr(c) = carr(c) + parr(p) * veg_pp%wtcol(p)
             sumwt(c) = sumwt(c) + veg_pp%wtcol(p)
          end if
       end if
    end do
    found = .false.
    do c = bounds%begc,bounds%endc
       if (sumwt(c) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = c
       else if (sumwt(c) /= 0._r8) then
          carr(c) = carr(c)/sumwt(c)
       end if
    end do
    if (found) then
      stop
    end if

  end subroutine p2c_1d_gpu

  subroutine p2c_1d_parallel (bounds, parr, carr, p2c_scale_type, para)
   ! !DESCRIPTION:
   ! Perfrom subgrid-average from pfts to columns.
   ! Averaging is only done for points that are not equal to "spval".
   !
   ! !ARGUMENTS:
   type(bounds_type), intent(in) :: bounds
   real(r8), intent(in)  :: parr( bounds%begp: )   ! patch array
   real(r8), intent(out) :: carr(1:)   ! column array
   integer, intent(in) :: p2c_scale_type ! scale type
   logical, intent(in) :: para 
   !
   ! !LOCAL VARIABLES:
   integer  :: p,c,index                       ! indices
   logical  :: found                              ! temporary for error check
   real(r8) :: sumwt, sum_arr !(bounds%begc:bounds%endc)     ! sum of weights
   !------------------------------------------------------------------------

   !$acc parallel loop independent gang worker default(present) private(sumwt,sum_arr)
   do c = bounds%begc, bounds%endc
      sumwt = 0._r8
      sum_arr = 0._r8
      !$acc loop vector reduction(+:sumwt)
      do p = col_pp%pfti(c), col_pp%pftf(c)
         if (veg_pp%active(p) .and. veg_pp%wtcol(p) /= 0._r8) then
            if (parr(p) /= spval) then
               
               sum_arr = sum_arr + parr(p) * veg_pp%wtcol(p)
               sumwt = sumwt + veg_pp%wtcol(p)
            end if
         end if
      end do
      carr(c) = sum_arr ! is this an issue for not being spval anymore?
      if(sumwt > 1.0_r8 + 1.e-6_r8 ) then 
         stop "Error p2c_1d sumwts > 1"
      elseif(sumwt /= 0._r8) then 
         carr(c) = carr(c)/sumwt
      elseif(sumwt == 0._r8) then 
         carr(c) = spval 
      end if 
   end do 
   
   ! found = .false.
   ! do c = bounds%begc,bounds%endc
   !    if (sumwt(c) > 1.0_r8 + 1.e-6_r8) then
   !       found = .true.
   !       index = c
   !    else if (sumwt(c) /= 0._r8) then
   !       carr(c) = carr(c)/sumwt(c)
   !    end if
   ! end do
   ! if (found) then
   !   stop
   ! end if
 end subroutine p2c_1d_parallel 

  !-----------------------------------------------------------------------
  subroutine p2c_2d (bounds, num2d, parr, carr, p2c_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from landunits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds
    integer           , intent(in)  :: num2d                     ! size of second dimension
    real(r8)          , intent(in)  :: parr( bounds%begp: , 1: ) ! patch array
    real(r8)          , intent(out) :: carr( bounds%begc: , 1: ) ! column array
    character(len=*)  , intent(in)  :: p2c_scale_type     ! scale type
    !
    ! !LOCAL VARIABLES:
    integer  :: j,p,c,index                         ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: sumwt(bounds%begc:bounds%endc)         ! sum of weights
    !------------------------------------------------------------------------

    carr(bounds%begc : bounds%endc, :) = spval
    do j = 1,num2d
       sumwt(bounds%begc : bounds%endc) = 0._r8
       do p = bounds%begp,bounds%endp
          if (veg_pp%active(p) .and. veg_pp%wtcol(p) /= 0._r8) then
             if (parr(p,j) /= spval) then
                c = veg_pp%column(p)
                if (sumwt(c) == 0._r8) carr(c,j) = 0._r8
                carr(c,j) = carr(c,j) + parr(p,j) * veg_pp%wtcol(p)
                sumwt(c) = sumwt(c) + veg_pp%wtcol(p)
             end if
          end if
       end do
       found = .false.
       do c = bounds%begc,bounds%endc
          if (sumwt(c) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = c
          else if (sumwt(c) /= 0._r8) then
             carr(c,j) = carr(c,j)/sumwt(c)
          end if
       end do

    end do
  end subroutine p2c_2d


  !-----------------------------------------------------------------------
  subroutine p2c_2d_parallel(bounds, num2d, parr, carr)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from landunits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds
    integer           , intent(in)  :: num2d                     ! size of second dimension
    real(r8)          , intent(in)  :: parr( bounds%begp: , 1: ) ! patch array
    real(r8)          , intent(out) :: carr( bounds%begc: , 1: ) ! column array
    !
    ! !LOCAL VARIABLES:
    integer  :: j,p,c,index                      ! indices
    real(r8), parameter :: scale_p2c = 1.0_r8    ! scale factor for column->landunit mapping
    logical  :: found                            ! temporary for error check
    real(r8) :: sumwt(bounds%begc:bounds%endc)   ! sum of weights
    real(r8) :: sum1 
    !------------------------------------------------------------------------

    !$acc enter data create(&
    !$acc sumwt(:), &
    !$acc sum1)

    !$acc parallel loop independent gang worker default(present) private(sum1)
    do c = bounds%begc, bounds%endc 
      sum1 = 0._r8 
      !$acc loop vector reduction(+:sum1)
      do p = col_pp%pfti(c), col_pp%pftf(c) 
         if (veg_pp%active(p) .and. veg_pp%wtcol(p) /= 0._r8 .and. parr(p,j) /= spval) then
            sum1 = sum1 + veg_pp%wtcol(p)
         end if 
      end do
      if(sum1 > 1.0_r8 + 1.e-6_r8) stop  
      sumwt(c) = sum1 
    end do

    ! carr(bounds%begc : bounds%endc, :) = spval
    !$acc parallel loop independent gang worker default(present) private(sum1)
    do j = 1,num2d
      do c = bounds%begc,bounds%endc
            sum1 = 0._r8 
         do p = col_pp%pfti(c), col_pp%pftf(c) 
            if (veg_pp%active(p) .and. veg_pp%wtcol(p) /= 0._r8 .and. parr(p,j) /= spval) then
                sum1 = sum1 + parr(p,j) * veg_pp%wtcol(p)
            end if
         end do
         carr(c,j) = sum1
         if(sumwt(c) /= 0._r8) then 
            carr(c,j) = carr(c,j)/sumwt(c)
         end if 
      end do 
    end do 

    !$acc exit data delete(&
    !$acc sumwt(:), &
    !$acc sum1)

  end subroutine p2c_2d_parallel


  !-----------------------------------------------------------------------
  subroutine p2c_1d_filter_parallel( numfc, filterc,  pftarr, colarr)
    !
    ! !DESCRIPTION:
    ! perform pft to column averaging for single level pft arrays
    ! Note: pftarr and colarr are the entire array for a processor
    !      not divisble by clumps
    !
    ! !ARGUMENTS:
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
    real(r8), intent(in)  :: pftarr(:)
    real(r8), intent(out) :: colarr(:)

    ! !LOCAL VARIABLES:
    integer :: fc,c,p  ! indices
    real(r8) :: sum1
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes

    !$acc parallel loop independent gang worker private(c,sum1) &
    !$acc default(present) create(sum1)
    do fc = 1,numfc
       c = filterc(fc)
       sum1 = 0._r8
       !$acc loop vector reduction(+:sum1)
       do p = col_pp%pfti(c), col_pp%pftf(c)
          if (veg_pp%active(p)) sum1 = sum1 + pftarr(p) * veg_pp%wtcol(p)
       end do
       colarr(c) = sum1
    end do

end subroutine p2c_1d_filter_parallel

  !-----------------------------------------------------------------------
  subroutine p2c_1d_filter (bounds, numfc, filterc,  pftarr, colarr)
    !
    ! !DESCRIPTION:
    ! perform pft to column averaging for single level pft arrays
    !
    !$acc routine seq
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
    real(r8), intent(in)  :: pftarr( bounds%begp: )
    real(r8), intent(out) :: colarr( bounds%begc: )

       ! !LOCAL VARIABLES:
    integer :: fc,c,p  ! indices
    !-----------------------------------------------------------------------
    do fc = 1,numfc
       c = filterc(fc)
       colarr(c) = 0._r8
       do p = col_pp%pfti(c), col_pp%pftf(c)
          if (veg_pp%active(p)) colarr(c) = colarr(c) + pftarr(p) * veg_pp%wtcol(p)
       end do
    end do

  end subroutine p2c_1d_filter

  !-----------------------------------------------------------------------
  subroutine p2c_2d_filter (lev, numfc, filterc, pftarr, colarr)
    !$acc routine seq
    ! !DESCRIPTION:
    ! perform pft to column averaging for multi level pft arrays
    !
    ! !ARGUMENTS:
    integer , intent(in)  :: lev
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
    real(r8), pointer     :: pftarr(:,:)
    real(r8), pointer     :: colarr(:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,p,j    ! indices
    !-----------------------------------------------------------------------

    do j = 1,lev
       do fc = 1,numfc
          c = filterc(fc)
          colarr(c,j) = 0._r8
          do p = col_pp%pfti(c), col_pp%pftf(c)
             if (veg_pp%active(p)) colarr(c,j) = colarr(c,j) + pftarr(p,j) * veg_pp%wtcol(p)
          end do
       end do
    end do

  end subroutine p2c_2d_filter

  !-----------------------------------------------------------------------
  subroutine p2l_1d (bounds, parr, larr, p2c_scale_type, c2l_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to landunits
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )  ! input column array
    real(r8), intent(out) :: larr( bounds%begl: )  ! output landunit array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,l,index                     ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: sumwt(bounds%begl:bounds%endl)     ! sum of weights
    real(r8) :: scale_p2c(bounds%begc:bounds%endc) ! scale factor for pft->column mapping
    real(r8) :: scale_c2l(bounds%begc:bounds%endc) ! scale factor for column->landunit mapping
    !------------------------------------------------------------------------
    ! Enforce expected array sizes

    call create_scale_c2l(bounds,c2l_scale_type,scale_c2l(bounds%begc:bounds%endc))

    if (p2c_scale_type == 'unity') then
       do p = bounds%begp,bounds%endp
          scale_p2c(p) = 1.0_r8
       end do
    else
       !#py write(iulog,*)'p2l_1d error: scale type ',p2c_scale_type,' not supported'
       !#py !#py call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    larr(bounds%begl : bounds%endl) = spval
    sumwt(bounds%begl : bounds%endl) = 0._r8
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .and. veg_pp%wtlunit(p) /= 0._r8) then
          c = veg_pp%column(p)
          if (parr(p) /= spval .and. scale_c2l(c) /= spval) then
             l = veg_pp%landunit(p)
             if (sumwt(l) == 0._r8) larr(l) = 0._r8
             larr(l) = larr(l) + parr(p) * scale_p2c(p) * scale_c2l(c) * veg_pp%wtlunit(p)
             sumwt(l) = sumwt(l) + veg_pp%wtlunit(p)
          end if
       end if
    end do
    found = .false.
    do l = bounds%begl,bounds%endl
       if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = l
       else if (sumwt(l) /= 0._r8) then
          larr(l) = larr(l)/sumwt(l)
       end if
    end do
    if (found) then
       !#py write(iulog,*)'p2l_1d error: sumwt is greater than 1.0 at l= ',index
       !#py !#py call endrun(decomp_index=index, elmlevel=namel, msg=errMsg(__FILE__, __LINE__))
    end if

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
    integer , intent(in)  :: num2d                 ! size of second dimension
    real(r8), intent(in)  :: parr( bounds%begp: , 1: )  ! input pft array
    real(r8), intent(out) :: larr( bounds%begl: , 1: )  ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,p,c,l,index       ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: sumwt(bounds%begl:bounds%endl)         ! sum of weights
    real(r8) :: scale_p2c(bounds%begc:bounds%endc)     ! scale factor for pft->column mapping
    real(r8) :: scale_c2l(bounds%begc:bounds%endc)     ! scale factor for column->landunit mapping
    !------------------------------------------------------------------------

    ! Enforce expected array sizes
    
    call create_scale_c2l(bounds,c2l_scale_type,scale_c2l(bounds%begc:bounds%endc))


    if (p2c_scale_type == 'unity') then
       do p = bounds%begp,bounds%endp
          scale_p2c(p) = 1.0_r8
       end do
    else
       !#py write(iulog,*)'p2l_2d error: scale type ',p2c_scale_type,' not supported'
       !#py !#py call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    larr(bounds%begl : bounds%endl, :) = spval
    do j = 1,num2d
       sumwt(bounds%begl : bounds%endl) = 0._r8
       do p = bounds%begp,bounds%endp
          if (veg_pp%active(p) .and. veg_pp%wtlunit(p) /= 0._r8) then
             c = veg_pp%column(p)
             if (parr(p,j) /= spval .and. scale_c2l(c) /= spval) then
                l = veg_pp%landunit(p)
                if (sumwt(l) == 0._r8) larr(l,j) = 0._r8
                larr(l,j) = larr(l,j) + parr(p,j) * scale_p2c(p) * scale_c2l(c) * veg_pp%wtlunit(p)
                sumwt(l) = sumwt(l) + veg_pp%wtlunit(p)
             end if
          end if
       end do
       found = .false.
       do l = bounds%begl,bounds%endl
          if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = l
          else if (sumwt(l) /= 0._r8) then
             larr(l,j) = larr(l,j)/sumwt(l)
          end if
       end do
       if (found) then
          !#py write(iulog,*)'p2l_2d error: sumwt is greater than 1.0 at l= ',index,' j= ',j
          !#py !#py call endrun(decomp_index=index, elmlevel=namel, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine p2l_2d

  !-----------------------------------------------------------------------
  subroutine p2g_1d(bounds, parr, garr, p2c_scale_type, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )  ! input pft array
    real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
    !
    !  !LOCAL VARIABLES:
    integer  :: p,c,l,g,index                   ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_p2c(bounds%begp:bounds%endp) ! scale factor
    real(r8) :: scale_c2l(bounds%begc:bounds%endc) ! scale factor
    real(r8) :: scale_l2g(bounds%begl:bounds%endl) ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    call build_scale_l2g(bounds, l2g_scale_type, &
         scale_l2g)

    call create_scale_c2l(bounds,c2l_scale_type,scale_c2l(bounds%begc:bounds%endc))

    if (p2c_scale_type == 'unity') then
       do p = bounds%begp,bounds%endp
          scale_p2c(p) = 1.0_r8
       end do
    else
       !#py write(iulog,*)'p2g_1d error: scale type ',c2l_scale_type,' not supported'
       !#py !#py call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    garr(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .and. veg_pp%wtgcell(p) /= 0._r8) then
          c = veg_pp%column(p)
          l = veg_pp%landunit(p)
          if (parr(p) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
             g = veg_pp%gridcell(p)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + parr(p) * scale_p2c(p) * scale_c2l(c) * scale_l2g(l) * veg_pp%wtgcell(p)
             sumwt(g) = sumwt(g) + veg_pp%wtgcell(p)
          end if
       end if
    end do
    found = .false.
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       !#py write(iulog,*)'p2g_1d error: sumwt is greater than 1.0 at g= ',index
       !#py !#py call endrun(decomp_index=index, elmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine p2g_1d

  !-----------------------------------------------------------------------
  subroutine p2g_1d_gpu(bounds, parr, garr, p2c_scale_type, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !$acc routine seq
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: parr( bounds%begp: )  ! input pft array
    real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
    integer , intent(in) :: p2c_scale_type ! scale factor type for averaging
    integer , intent(in) :: c2l_scale_type ! scale factor type for averaging |||  unity = 0, urbanf = 1, urbans = 2
    integer , intent(in) :: l2g_scale_type ! unity =0, natveg = 3, veg =4, ice=5, nonurb=6, lake=7
    !
    !  !LOCAL VARIABLES:
    integer  :: p,c,l,g,index                   ! indices
    logical  :: found                              ! temporary for error check
    real(r8),parameter :: scale_p2c = 1._r8 ! scale factor
    real(r8) :: scale_c2l ! scale factor
    real(r8) :: scale_l2g ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------
    garr(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .and. veg_pp%wtgcell(p) /= 0._r8) then
          c = veg_pp%column(p)
          l = veg_pp%landunit(p)
          scale_l2g = main_scale_l2g_lookup(l2g_scale_type,lun_pp%itype(l))
          scale_c2l = main_scale_c2l(c,c2l_scale_type)
          if (parr(p) /= spval .and. scale_c2l /= spval .and. scale_l2g /= spval) then
             g = veg_pp%gridcell(p)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + parr(p) * scale_p2c * scale_c2l * scale_l2g * veg_pp%wtgcell(p)
             sumwt(g) = sumwt(g) + veg_pp%wtgcell(p)
          end if
       end if
    end do
    found = .false.
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       print *,'p2g_1d error: sumwt is greater than 1.0 at g= ',index
    end if

  end subroutine p2g_1d_gpu


  !-----------------------------------------------------------------------
  subroutine p2g_2d(bounds, num2d, parr, garr, p2c_scale_type, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: parr( bounds%begp: , 1: ) ! input pft array
    real(r8), intent(out) :: garr( bounds%begg: , 1: ) ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type     ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type     ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type     ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,p,c,l,g,index                     ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_p2c(bounds%begp:bounds%endp)     ! scale factor
    real(r8) :: scale_c2l(bounds%begc:bounds%endc)     ! scale factor
    real(r8) :: scale_l2g(bounds%begl:bounds%endl)     ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)         ! sum of weights
    !------------------------------------------------------------------------
    ! Enforce expected array sizes

    call build_scale_l2g(bounds, l2g_scale_type, &
         scale_l2g)

    call create_scale_c2l(bounds,c2l_scale_type, scale_c2l(bounds%begc:bounds%endc))

    if (p2c_scale_type == 'unity') then
       do p = bounds%begp,bounds%endp
          scale_p2c(p) = 1.0_r8
       end do
    else
       !#py write(iulog,*)'p2g_2d error: scale type ',c2l_scale_type,' not supported'
       !#py !#py call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    garr(bounds%begg : bounds%endg, :) = spval
    do j = 1,num2d
       sumwt(bounds%begg : bounds%endg) = 0._r8
       do p = bounds%begp,bounds%endp
          if (veg_pp%active(p) .and. veg_pp%wtgcell(p) /= 0._r8) then
             c = veg_pp%column(p)
             l = veg_pp%landunit(p)
             if (parr(p,j) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
                g = veg_pp%gridcell(p)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + parr(p,j) * scale_p2c(p) * scale_c2l(c) * scale_l2g(l) * veg_pp%wtgcell(p)
                sumwt(g) = sumwt(g) + veg_pp%wtgcell(p)
             end if
          end if
       end do
       found = .false.
       do g = bounds%begg, bounds%endg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          !#py write(iulog,*)'p2g_2d error: sumwt gt 1.0 at g/sumwt = ',index,sumwt(index)
          !#py !#py call endrun(decomp_index=index, elmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine p2g_2d

  !-----------------------------------------------------------------------
  subroutine p2g_2d_gpu(bounds, num2d, parr, garr, p2c_scale_type, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    !$acc routine seq
    ! Perfrom subgrid-average from pfts to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: parr( bounds%begp: , 1: ) ! input pft array
    real(r8), intent(out) :: garr( bounds%begg: , 1: ) ! output gridcell array
    integer, intent(in) :: p2c_scale_type     ! scale factor type for averaging
    integer, intent(in) :: c2l_scale_type     ! scale factor type for averaging
    integer, intent(in) :: l2g_scale_type     ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,p,c,l,g,index                     ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_p2c(bounds%begp:bounds%endp)     ! scale factor
    real(r8) :: scale_c2l(bounds%begc:bounds%endc)     ! scale factor
    real(r8) :: scale_l2g(bounds%begl:bounds%endl)     ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)         ! sum of weights
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    call build_scale_l2g_gpu(bounds, l2g_scale_type, &
         scale_l2g)

    if (c2l_scale_type == unity) then
       do c = bounds%begc,bounds%endc
          scale_c2l(c) = 1.0_r8
       end do
    else if (c2l_scale_type == urbanf) then
       do c = bounds%begc,bounds%endc
          l = col_pp%landunit(c)
          if (lun_pp%urbpoi(l)) then
             if (col_pp%itype(c) == icol_sunwall) then
                scale_c2l(c) = 3.0 * lun_pp%canyon_hwr(l)
             else if (col_pp%itype(c) == icol_shadewall) then
                scale_c2l(c) = 3.0 * lun_pp%canyon_hwr(l)
             else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0_r8
             else if (col_pp%itype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else if (c2l_scale_type == urbans) then
       do c = bounds%begc,bounds%endc
          l = col_pp%landunit(c)
          if (lun_pp%urbpoi(l)) then
             if (col_pp%itype(c) == icol_sunwall) then
                scale_c2l(c) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
             else if (col_pp%itype(c) == icol_shadewall) then
                scale_c2l(c) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
             else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0 / (2.*lun_pp%canyon_hwr(l) + 1.)
             else if (col_pp%itype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else
       !write(iulog,*)'p2g_2d error: scale type ',c2l_scale_type,' not supported'
       !call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (p2c_scale_type == unity) then
       do p = bounds%begp,bounds%endp
          scale_p2c(p) = 1.0_r8
       end do
    else
       !write(iulog,*)'p2g_2d error: scale type ',c2l_scale_type,' not supported'
       !call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    garr(bounds%begg : bounds%endg, :) = spval
    do j = 1,num2d
       sumwt(bounds%begg : bounds%endg) = 0._r8
       do p = bounds%begp,bounds%endp
          if (veg_pp%active(p) .and. veg_pp%wtgcell(p) /= 0._r8) then
             c = veg_pp%column(p)
             l = veg_pp%landunit(p)
             if (parr(p,j) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
                g = veg_pp%gridcell(p)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + parr(p,j) * scale_p2c(p) * scale_c2l(c) * scale_l2g(l) * veg_pp%wtgcell(p)
                sumwt(g) = sumwt(g) + veg_pp%wtgcell(p)
             end if
          end if
       end do
       found = .false.
       do g = bounds%begg, bounds%endg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          print *,'p2g_2d error: sumwt gt 1.0 at g/sumwt = ',index,sumwt(index)
          stop
       end if
    end do

  end subroutine p2g_2d_gpu

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
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,index                       ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_c2l(bounds%begc:bounds%endc) ! scale factor for column->landunit mapping
    real(r8) :: sumwt(bounds%begl:bounds%endl)     ! sum of weights
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    call create_scale_c2l(bounds, c2l_scale_type, scale_c2l(bounds%begc:bounds%endc))

    larr(bounds%begl : bounds%endl) = spval
    sumwt(bounds%begl : bounds%endl) = 0._r8
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .and. col_pp%wtlunit(c) /= 0._r8) then
          if (carr(c) /= spval .and. scale_c2l(c) /= spval) then
             l = col_pp%landunit(c)
             if (sumwt(l) == 0._r8) larr(l) = 0._r8
             larr(l) = larr(l) + carr(c) * scale_c2l(c) * col_pp%wtlunit(c)
             sumwt(l) = sumwt(l) + col_pp%wtlunit(c)
          end if
       end if
    end do
    found = .false.
    do l = bounds%begl,bounds%endl
       if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = l
       else if (sumwt(l) /= 0._r8) then
          larr(l) = larr(l)/sumwt(l)
       end if
    end do
    if (found) then
       !#py write(iulog,*)'c2l_1d error: sumwt is greater than 1.0 at l= ',index
       !#py !#py call endrun(decomp_index=index, elmlevel=namel, msg=errMsg(__FILE__, __LINE__))
    end if

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
    character(len=*), intent(in) :: c2l_scale_type     ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,index                         ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_c2l(bounds%begc:bounds%endc)     ! scale factor for column->landunit mapping
    real(r8) :: sumwt(bounds%begl:bounds%endl)         ! sum of weights
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    call create_scale_c2l(bounds, c2l_scale_type, &
         scale_c2l(bounds%begc:bounds%endc))
    
    larr(bounds%begl : bounds%endl, :) = spval
    do j = 1,num2d
       sumwt(bounds%begl : bounds%endl) = 0._r8
       do c = bounds%begc,bounds%endc
          if (col_pp%active(c) .and. col_pp%wtlunit(c) /= 0._r8) then
             if (carr(c,j) /= spval .and. scale_c2l(c) /= spval) then
                l = col_pp%landunit(c)
                if (sumwt(l) == 0._r8) larr(l,j) = 0._r8
                larr(l,j) = larr(l,j) + carr(c,j) * scale_c2l(c) * col_pp%wtlunit(c)
                sumwt(l) = sumwt(l) + col_pp%wtlunit(c)
             end if
          end if
       end do
       found = .false.
       do l = bounds%begl,bounds%endl
          if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = l
          else if (sumwt(l) /= 0._r8) then
             larr(l,j) = larr(l,j)/sumwt(l)
          end if
       end do
       if (found) then
        !#py write(iulog,*)'c2l_2d error: sumwt is greater than 1.0 at l= ',index,' lev= ',j
        !#py !#py call endrun(decomp_index=index, elmlevel=namel, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine c2l_2d

  !-----------------------------------------------------------------------
  subroutine c2g_1d(bounds, carr, garr, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from columns to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: carr( bounds%begc: )  ! input column array
    real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,g,index                     ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_c2l(bounds%begc:bounds%endc) ! scale factor
    real(r8) :: scale_l2g(bounds%begl:bounds%endl) ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------

    call build_scale_l2g(bounds, l2g_scale_type, &
         scale_l2g)

    call create_scale_c2l(bounds,c2l_scale_type, scale_c2l(bounds%begc:bounds%endc))

    garr(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .and. col_pp%wtgcell(c) /= 0._r8) then
          l = col_pp%landunit(c)
          if (carr(c) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
             g = col_pp%gridcell(c)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + carr(c) * scale_c2l(c) * scale_l2g(l) * col_pp%wtgcell(c)
             sumwt(g) = sumwt(g) + col_pp%wtgcell(c)
          end if
       end if
    end do
    found = .false.
    
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       !#py write(iulog,*)'c2g_1d error: sumwt is greater than 1.0 at g= ',index
       !#py !#py call endrun(decomp_index=index, elmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine c2g_1d
  !-----------------------------------------------------------------------
  subroutine c2g_1d_gpu(bounds, carr, garr, c2l_scale_type, l2g_scale_type)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from columns to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: carr( bounds%begc: )  ! input column array
    real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
    integer , intent(in) :: c2l_scale_type  !! unity = 0, urbanf = 1, urbans = 2
    integer , intent(in) :: l2g_scale_type  !!natveg = 3, veg =4, ice=5, nonurb=6, lake=7
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,g,index                     ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_c2l ! scale factor
    real(r8) :: scale_l2g ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------
   
    garr(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8

    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .and. col_pp%wtgcell(c) /= 0._r8) then
          l = col_pp%landunit(c)
          scale_l2g = main_scale_l2g_lookup(l2g_scale_type,lun_pp%itype(l))
          scale_c2l = main_scale_c2l(c,c2l_scale_type)
          if (carr(c) /= spval .and. scale_c2l /= spval .and. scale_l2g /= spval) then
            g = col_pp%gridcell(c)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + carr(c) * scale_c2l * scale_l2g * col_pp%wtgcell(c)
             sumwt(g) = sumwt(g) + col_pp%wtgcell(c)
          end if
       end if
    end do
    found = .false.
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do

  end subroutine c2g_1d_gpu

  
  subroutine c2g_1d_parallel(bounds, carr, garr, c2l_scale_type, l2g_scale_type,para)
   ! !DESCRIPTION:
   ! Perfrom subgrid-average from columns to gridcells.
   ! Averaging is only done for points that are not equal to "spval".
   !
   use GridcellType , only : grc_pp 
   ! !ARGUMENTS:
   type(bounds_type), intent(in) :: bounds
   real(r8), intent(in)  :: carr( bounds%begc: )  ! input column array
   real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
   integer , intent(in) :: c2l_scale_type  !! unity = 1, urbanf = 2, urbans = 3
   integer , intent(in) :: l2g_scale_type  !!natveg = 2, veg =3, ice=4, nonurb=5, lake=6
   logical, intent(in) :: para
   !
   ! !LOCAL VARIABLES:
   integer  :: c,l,g,index,fc  ! indices
   logical  :: found        ! temporary for error check
   real(r8) :: scale_c2l  ! scale factor ! now using main_scale_c2l
   real(r8) :: scale_l2g  ! scale factor
   real(r8) :: sumwt   ! sum of weights
   real(r8) :: sum_g 
   !- -----------------------------------------------------------------------
   ! note : scale_l2g(l) = main_scale_l2g_lookup(TYPE,lun_pp%itype(l))
   !        main_scale_c2l(cols, c2l_scale_type)
   !$acc enter data create(&
   !$acc sumwt, &
   !$acc sum_g)

   !$acc parallel loop independent gang worker default(present) private(sumwt,sum_g) copyin(c2l_scale_type)
   do g = bounds%begg,bounds%endg
      sumwt = 0._r8
      sum_g = 0._r8
      !$acc loop reduction(+:sumwt,sum_g) private(l,scale_l2g,scale_c2l,c)
      do fc = 1, grc_pp%ncolumns(g)
         c = grc_pp%cols(g,fc)  
         if (col_pp%active(c) .and. col_pp%wtgcell(c) /= 0._r8) then
            l = col_pp%landunit(c)
            scale_l2g = main_scale_l2g_lookup(l2g_scale_type,lun_pp%itype(l))
            scale_c2l = main_scale_c2l(c,c2l_scale_type)
            if (carr(c) /= spval .and. scale_c2l /= spval .and. scale_l2g /= spval) then
               sum_g = sum_g + carr(c) * scale_c2l * scale_l2g * col_pp%wtgcell(c)
               sumwt = sumwt + col_pp%wtgcell(c)
            end if
         end if
      end do 
      garr(g) = sum_g

      if(sumwt .ne. 0._r8) then 
         garr(g) = garr(g)/sumwt
      elseif(sumwt == 0._r8 ) then 
         garr(g) = spval  !! needed to keep unused equal to spval?
      elseif(sumwt > 1.0_r8 + 1.e-6_r8 ) then 
         stop "Error with col gridcell weights"
      end if 
   end do

   !$acc exit data delete(&
   !$acc sumwt, &
   !$acc sum_g)

 end subroutine c2g_1d_parallel

  !-----------------------------------------------------------------------
  subroutine c2g_2d(bounds, num2d, carr, garr, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from columns to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: carr( bounds%begc: , 1: ) ! input column array
    real(r8), intent(out) :: garr( bounds%begg: , 1: ) ! output gridcell array
    character(len=*), intent(in) :: c2l_scale_type     ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type     ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,g,l,index                       ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_c2l(bounds%begc:bounds%endc)     ! scale factor
    real(r8) :: scale_l2g(bounds%begl:bounds%endl)     ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)         ! sum of weights
    !------------------------------------------------------------------------

    call build_scale_l2g(bounds, l2g_scale_type, &
         scale_l2g)

    call create_scale_c2l(bounds, c2l_scale_type, scale_c2l(bounds%begc:bounds%endc))

    garr(bounds%begg : bounds%endg,:) = spval
    do j = 1,num2d
       sumwt(bounds%begg : bounds%endg) = 0._r8
       do c = bounds%begc,bounds%endc
          if (col_pp%active(c) .and. col_pp%wtgcell(c) /= 0._r8) then
             l = col_pp%landunit(c)
             if (carr(c,j) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
                g = col_pp%gridcell(c)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + carr(c,j) * scale_c2l(c) * scale_l2g(l) * col_pp%wtgcell(c)
                sumwt(g) = sumwt(g) + col_pp%wtgcell(c)
             end if
          end if
       end do
       found = .false.
       do g = bounds%begg, bounds%endg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          !#py write(iulog,*)'c2g_2d error: sumwt is greater than 1.0 at g= ',index
          !#py !#py call endrun(decomp_index=index, elmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine c2g_2d

  !-----------------------------------------------------------------------
  subroutine c2g_2d_gpu(bounds, num2d, carr, garr, c2l_scale_type, l2g_scale_type)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from columns to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: carr( bounds%begc: , 1: ) ! input column array
    real(r8), intent(out) :: garr( bounds%begg: , 1: ) ! output gridcell array
    integer , intent(in)  :: c2l_scale_type
    integer , intent(in)  :: l2g_scale_type ! unity =0, natveg = 3, veg =4, ice=5, nonurb=6, lake=7
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,g,l,index                       ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_c2l(bounds%begc:bounds%endc)     ! scale factor
    real(r8) :: scale_l2g(bounds%begl:bounds%endl)     ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)         ! sum of weights
    !------------------------------------------------------------------------

    call build_scale_l2g_gpu(bounds, l2g_scale_type, &
         scale_l2g)

   call create_scale_c2l_gpu(bounds,c2l_scale_type, scale_c2l)

    garr(bounds%begg : bounds%endg,:) = spval
    do j = 1,num2d
       sumwt(bounds%begg : bounds%endg) = 0._r8
       do c = bounds%begc,bounds%endc
          if (col_pp%active(c) .and. col_pp%wtgcell(c) /= 0._r8) then
             l = col_pp%landunit(c)
             if (carr(c,j) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
                g = col_pp%gridcell(c)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + carr(c,j) * scale_c2l(c) * scale_l2g(l) * col_pp%wtgcell(c)
                sumwt(g) = sumwt(g) + col_pp%wtgcell(c)
             end if
          end if
       end do
       found = .false.
       do g = bounds%begg, bounds%endg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          print *, 'c2g_2d error: sumwt is greater than 1.0 at g= ',index
       end if
    end do

  end subroutine c2g_2d_gpu

  !-----------------------------------------------------------------------
  subroutine l2g_1d(bounds, larr, garr, l2g_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from landunits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: larr( bounds%begl: )  ! input landunit array
    real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: l,g,index                       ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_l2g(bounds%begl:bounds%endl) ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------
    call build_scale_l2g(bounds, l2g_scale_type, &
         scale_l2g(bounds%begl:bounds%endl))

    garr(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8
    do l = bounds%begl,bounds%endl
       if (lun_pp%active(l) .and. lun_pp%wtgcell(l) /= 0._r8) then
          if (larr(l) /= spval .and. scale_l2g(l) /= spval) then
             g = lun_pp%gridcell(l)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + larr(l) * scale_l2g(l) * lun_pp%wtgcell(l)
             sumwt(g) = sumwt(g) + lun_pp%wtgcell(l)
          end if
       end if
    end do
    found = .false.
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       !#py write(iulog,*) 'l2g_1d error: sumwt is greater than 1.0 at g= ',index
       !#py !#py call endrun(decomp_index=index, elmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
    end if

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
    real(r8), intent(in)  :: larr( bounds%begl: , 1: ) ! input landunit array
    real(r8), intent(out) :: garr( bounds%begg: , 1: ) ! output gridcell array
    character(len=*), intent(in) :: l2g_scale_type     ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,g,l,index                         ! indices
    integer  :: max_lu_per_gcell                       ! max landunits per gridcell; on the fly
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_l2g(bounds%begl:bounds%endl)     ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)         ! sum of weights
    !------------------------------------------------------------------------

    call build_scale_l2g(bounds, l2g_scale_type, &
         scale_l2g(bounds%begl:bounds%endl))

    garr(bounds%begg : bounds%endg, :) = spval
    do j = 1,num2d
       sumwt(bounds%begg : bounds%endg) = 0._r8
       do l = bounds%begl,bounds%endl
          if (lun_pp%active(l) .and. lun_pp%wtgcell(l) /= 0._r8) then
             if (larr(l,j) /= spval .and. scale_l2g(l) /= spval) then
                g = lun_pp%gridcell(l)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + larr(l,j) * scale_l2g(l) * lun_pp%wtgcell(l)
                sumwt(g) = sumwt(g) + lun_pp%wtgcell(l)
             end if
          end if
       end do
       found = .false.
       do g = bounds%begg,bounds%endg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index= g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          !#py write(iulog,*) 'l2g_2d error: sumwt is greater than 1.0 at g= ',index,' lev= ',j
          !#py !#py call endrun(decomp_index=index, elmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine l2g_2d


  !-----------------------------------------------------------------------
  subroutine l2g_1d_gpu(bounds, larr, garr, l2g_scale_type)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from landunits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: larr( bounds%begl: )  ! input landunit array
    real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
    integer, intent(in) :: l2g_scale_type ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: l,g,index                       ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_l2g(bounds%begl:bounds%endl) ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------
    call build_scale_l2g_gpu(bounds, l2g_scale_type, &
         scale_l2g(bounds%begl:bounds%endl))

    garr(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8
    do l = bounds%begl,bounds%endl
       if (lun_pp%active(l) .and. lun_pp%wtgcell(l) /= 0._r8) then
          if (larr(l) /= spval .and. scale_l2g(l) /= spval) then
             g = lun_pp%gridcell(l)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + larr(l) * scale_l2g(l) * lun_pp%wtgcell(l)
             sumwt(g) = sumwt(g) + lun_pp%wtgcell(l)
          end if
       end if
    end do
    found = .false.
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       print *, 'l2g_1d error: sumwt is greater than 1.0 at g= ',index
    end if

  end subroutine l2g_1d_gpu

  !-----------------------------------------------------------------------
  subroutine l2g_2d_gpu(bounds, num2d, larr, garr, l2g_scale_type)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from landunits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: larr( bounds%begl: , 1: ) ! input landunit array
    real(r8), intent(out) :: garr( bounds%begg: , 1: ) ! output gridcell array
    integer, intent(in) :: l2g_scale_type     ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,g,l,index                         ! indices
    integer  :: max_lu_per_gcell                       ! max landunits per gridcell; on the fly
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_l2g(bounds%begl:bounds%endl)     ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)         ! sum of weights
    !------------------------------------------------------------------------

    call build_scale_l2g_gpu(bounds, l2g_scale_type, &
         scale_l2g(bounds%begl:bounds%endl))

    garr(bounds%begg : bounds%endg, :) = spval
    do j = 1,num2d
       sumwt(bounds%begg : bounds%endg) = 0._r8
       do l = bounds%begl,bounds%endl
          if (lun_pp%active(l) .and. lun_pp%wtgcell(l) /= 0._r8) then
             if (larr(l,j) /= spval .and. scale_l2g(l) /= spval) then
                g = lun_pp%gridcell(l)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + larr(l,j) * scale_l2g(l) * lun_pp%wtgcell(l)
                sumwt(g) = sumwt(g) + lun_pp%wtgcell(l)
             end if
          end if
       end do
       found = .false.
       do g = bounds%begg,bounds%endg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index= g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
         stop
       end if
    end do

  end subroutine l2g_2d_gpu
  
  
  !-----------------------------------------------------------------------
  subroutine build_scale_l2g(bounds, l2g_scale_type, scale_l2g)
    ! !DESCRIPTION:
    ! Fill the scale_l2g(bounds%begl:bounds%endl) array with appropriate values for the given l2g_scale_type.
    ! This array can later be used to scale each landunit in forming grid cell averages.
    !
    ! !USES:
    use landunit_varcon, only : max_lunit
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    character(len=*), intent(in)  :: l2g_scale_type            ! scale factor type for averaging
    real(r8)        , intent(out) :: scale_l2g( bounds%begl: ) ! scale factor
    !
    ! !LOCAL VARIABLES:
    real(r8) :: scale_lookup(max_lunit) ! scale factor for each landunit type
    integer  :: l                       ! index
    !-----------------------------------------------------------------------

     call create_scale_l2g_lookup(l2g_scale_type, scale_lookup)

     do l = bounds%begl,bounds%endl
        scale_l2g(l) = scale_lookup(lun_pp%itype(l))
     end do

  end subroutine build_scale_l2g

  !-----------------------------------------------------------------------
  subroutine build_scale_l2g_gpu(bounds, l2g_scale_type, scale_l2g)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Fill the scale_l2g(bounds%begl:bounds%endl) array with appropriate values for the given l2g_scale_type.
    ! This array can later be used to scale each landunit in forming grid cell averages.
    !
    ! !USES:
    use landunit_varcon, only : max_lunit
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !character(len=*), intent(in)  :: l2g_scale_type            ! scale factor type for averaging
    integer , intent(in)    :: l2g_scale_type ! unity =0, natveg = 3, veg =4, ice=5, nonurb=6, lake=7
    real(r8)        , intent(out) :: scale_l2g( bounds%begl: ) ! scale factor
    !
    ! !LOCAL VARIABLES:
    real(r8) :: scale_lookup(max_lunit) ! scale factor for each landunit type
    integer  :: l                       ! index
    !-----------------------------------------------------------------------


     call create_scale_l2g_lookup_gpu(l2g_scale_type, scale_lookup)

     do l = bounds%begl,bounds%endl
        scale_l2g(l) = scale_lookup(lun_pp%itype(l))
     end do

  end subroutine build_scale_l2g_gpu

  !-----------------------------------------------------------------------
  subroutine build_scale_l2t(bounds, l2t_scale_type, scale_l2t)
    !
    ! !DESCRIPTION:
    ! Fill the scale_l2t(bounds%begl:bounds%endl) array with appropriate values for the given l2t_scale_type for the topounit based structure
    ! This array can later be used to scale each landunit in forming grid cell averages.
    !
    ! !USES:
    use landunit_varcon, only : max_lunit
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds                    
    character(len=*), intent(in)  :: l2t_scale_type            ! scale factor type for averaging
    real(r8)        , intent(out) :: scale_l2t( bounds%begl: ) ! scale factor 
    !
    ! !LOCAL VARIABLES:
    real(r8) :: scale_lookup(max_lunit) ! scale factor for each landunit type
    integer  :: l                       ! index
    !-----------------------------------------------------------------------
     

     call create_scale_l2t_lookup(l2t_scale_type, scale_lookup)

     do l = bounds%begl,bounds%endl
        scale_l2t(l) = scale_lookup(lun_pp%itype(l))
     end do

  end subroutine build_scale_l2t
  
  !-----------------------------------------------------------------------
  subroutine create_scale_l2g_lookup(l2g_scale_type, scale_lookup)
    ! DESCRIPTION:
    ! Create a lookup array, scale_lookup(1..max_lunit), which gives the scale factor for
    ! each landunit type depending on l2g_scale_type
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop, istice, istice_mec, istdlak
    use landunit_varcon, only : isturb_MIN, isturb_MAX, max_lunit
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)  :: l2g_scale_type           ! scale factor type for averaging
    real(r8)   , intent(out) :: scale_lookup(max_lunit)  ! scale factor for each landunit type
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


     ! Initialize scale_lookup to spval for all landunits. Thus, any landunit that keeps
     ! the default value will be excluded from grid cell averages.
     scale_lookup(:) = spval

     if (trim(l2g_scale_type) == 'unity') then
        scale_lookup(:) = 1.0_r8
     else if (trim(l2g_scale_type) == 'natveg') then
        scale_lookup(istsoil) = 1.0_r8
     else if (trim(l2g_scale_type) == 'veg') then
        scale_lookup(istsoil) = 1.0_r8
        scale_lookup(istcrop) = 1.0_r8
     else if (trim(l2g_scale_type) == 'ice') then
        scale_lookup(istice) = 1.0_r8
        scale_lookup(istice_mec) = 1.0_r8
     else if (trim(l2g_scale_type) == 'nonurb') then
        scale_lookup(:) = 1.0_r8
        scale_lookup(isturb_MIN:isturb_MAX) = spval
     else if (trim(l2g_scale_type) == 'lake') then
        scale_lookup(istdlak) = 1.0_r8
     else
      !   write(iulog,*)'scale_l2g_lookup_array error: scale type ',l2g_scale_type,' not supported'
      !   call endrun(msg=errMsg(__FILE__, __LINE__))
     end if

  end subroutine create_scale_l2g_lookup
  
  !-----------------------------------------------------------------------
  subroutine create_scale_l2t_lookup(l2t_scale_type, scale_lookup)
    ! 
    ! DESCRIPTION:
    ! Create a lookup array, scale_lookup(1..max_lunit), which gives the scale factor for
    ! each landunit type depending on l2t_scale_type for the topounit based structure the same way to l2g_scale_type TKT
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop, istice, istice_mec, istdlak
    use landunit_varcon, only : isturb_MIN, isturb_MAX, max_lunit
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)  :: l2t_scale_type           ! scale factor type for averaging
    real(r8)        , intent(out) :: scale_lookup(max_lunit)  ! scale factor for each landunit type
    !-----------------------------------------------------------------------

     ! Initialize scale_lookup to spval for all landunits. Thus, any landunit that keeps
     ! the default value will be excluded from grid cell averages.
     scale_lookup(:) = spval

     if (l2t_scale_type == 'unity') then
        scale_lookup(:) = 1.0_r8
     else if (l2t_scale_type == 'natveg') then
        scale_lookup(istsoil) = 1.0_r8
     else if (l2t_scale_type == 'veg') then
        scale_lookup(istsoil) = 1.0_r8
        scale_lookup(istcrop) = 1.0_r8
     else if (l2t_scale_type == 'ice') then
        scale_lookup(istice) = 1.0_r8
        scale_lookup(istice_mec) = 1.0_r8
     else if (l2t_scale_type == 'nonurb') then
        scale_lookup(:) = 1.0_r8
        scale_lookup(isturb_MIN:isturb_MAX) = spval
     else if (l2t_scale_type == 'lake') then
        scale_lookup(istdlak) = 1.0_r8
     else
        !#py write(iulog,*)'scale_l2g_lookup_array error: scale type ',l2t_scale_type,' not supported'
        !#py !#py call endrun(msg=errMsg(__FILE__, __LINE__))
     end if

  end subroutine create_scale_l2t_lookup

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
    real(r8), intent(out) :: tarr( bounds%begt: )  ! output topounits array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2t_scale_type ! scale factor type for averaging
    !
    !  !LOCAL VARIABLES:
    integer  :: p,c,l,t,index                   ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_p2c(bounds%begp:bounds%endp) ! scale factor
    real(r8) :: scale_c2l(bounds%begc:bounds%endc) ! scale factor
    real(r8) :: scale_l2t(bounds%begl:bounds%endl) ! scale factor
    real(r8) :: sumwt(bounds%begt:bounds%endt)     ! sum of weights
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    call build_scale_l2t(bounds, l2t_scale_type, &
         scale_l2t(bounds%begl:bounds%endl))

    ! Build scale_c2l: TekluTesfa@PNNL created a new subroutine to remove repretitive code
    call create_scale_c2l(bounds, c2l_scale_type, &
         scale_c2l(bounds%begc:bounds%endc))
    
    if (p2c_scale_type == 'unity') then
       do p = bounds%begp,bounds%endp
          scale_p2c(p) = 1.0_r8
       end do
    else
       !#py write(iulog,*)'p2t_1d error: scale type ',p2c_scale_type,' not supported'
       !#py !#py call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    tarr(bounds%begt : bounds%endt) = spval
    sumwt(bounds%begt : bounds%endt) = 0._r8
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .and. veg_pp%wttopounit(p) /= 0._r8) then
          c = veg_pp%column(p)
          l = veg_pp%landunit(p)
          if (parr(p) /= spval .and. scale_c2l(c) /= spval .and. scale_l2t(l) /= spval) then
             t = veg_pp%topounit(p)
             if (sumwt(t) == 0._r8) tarr(t) = 0._r8
             tarr(t) = tarr(t) + parr(p) * scale_p2c(p) * scale_c2l(c) * scale_l2t(l) * veg_pp%wttopounit(p)
             sumwt(t) = sumwt(t) + veg_pp%wttopounit(p)
          end if
       end if
    end do
    found = .false.
    do t = bounds%begt, bounds%endt
       if (sumwt(t) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = t
       else if (sumwt(t) /= 0._r8) then
          tarr(t) = tarr(t)/sumwt(t)
       end if
    end do
    if (found) then
       !#py write(iulog,*)'p2t_1d error: sumwt is greater than 1.0 at t= ',index
       !#py !#py call endrun(decomp_index=index, elmlevel=namet, msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine p2t_1d

  !-----------------------------------------------------------------------
  subroutine p2t_2d(bounds, num2d, parr, tarr, p2c_scale_type, c2l_scale_type, l2t_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds            
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: parr( bounds%begp: , 1: ) ! input pft array
    real(r8), intent(out) :: tarr( bounds%begt: , 1: ) ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type     ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type     ! scale factor type for averaging
    character(len=*), intent(in) :: l2t_scale_type     ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,p,c,l,t,index                     ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_p2c(bounds%begp:bounds%endp)     ! scale factor
    real(r8) :: scale_c2l(bounds%begc:bounds%endc)     ! scale factor
    real(r8) :: scale_l2t(bounds%begl:bounds%endl)     ! scale factor
    real(r8) :: sumwt(bounds%begt:bounds%endt)         ! sum of weights
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    call build_scale_l2t(bounds, l2t_scale_type, &
         scale_l2t(bounds%begl:bounds%endl))

    ! Build scale_c2l: TekluTesfa@PNNL created a new subroutine to remove repretitive code
    call create_scale_c2l(bounds, c2l_scale_type, &
         scale_c2l(bounds%begc:bounds%endc))
    
    if (p2c_scale_type == 'unity') then
       do p = bounds%begp,bounds%endp
          scale_p2c(p) = 1.0_r8
       end do
    else
       !#py write(iulog,*)'p2t_2d error: scale type ',p2c_scale_type,' not supported'
       !#py !#py call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    tarr(bounds%begt : bounds%endt, :) = spval
    do j = 1,num2d
       sumwt(bounds%begt : bounds%endt) = 0._r8
       do p = bounds%begp,bounds%endp 
          if (veg_pp%active(p) .and. veg_pp%wttopounit(p) /= 0._r8) then
             c = veg_pp%column(p)
             l = veg_pp%landunit(p)
             if (parr(p,j) /= spval .and. scale_c2l(c) /= spval .and. scale_l2t(l) /= spval) then
                t = veg_pp%topounit(p)
                if (sumwt(t) == 0._r8) tarr(t,j) = 0._r8
                tarr(t,j) = tarr(t,j) + parr(p,j) * scale_p2c(p) * scale_c2l(c) * scale_l2t(l) * veg_pp%wttopounit(p)
                sumwt(t) = sumwt(t) + veg_pp%wttopounit(p)
             end if
          end if
       end do
       found = .false.
       do t = bounds%begt, bounds%endt
          if (sumwt(t) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = t
          else if (sumwt(t) /= 0._r8) then
             tarr(t,j) = tarr(t,j)/sumwt(t)
          end if
       end do
       if (found) then
          !#py write(iulog,*)'p2t_2d error: sumwt gt 1.0 at t/sumwt = ',index,sumwt(index)
          !#py !#py call endrun(decomp_index=index, elmlevel=namet, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

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
    real(r8), intent(out) :: tarr( bounds%begt: )  ! output gridcell array
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2t_scale_type ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,t,index                     ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_c2l(bounds%begc:bounds%endc) ! scale factor
    real(r8) :: scale_l2t(bounds%begl:bounds%endl) ! scale factor
    real(r8) :: sumwt(bounds%begt:bounds%endt)     ! sum of weights
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    call build_scale_l2t(bounds, l2t_scale_type, &
         scale_l2t(bounds%begl:bounds%endl))

    ! Build scale_c2l: TekluTesfa@PNNL created a new subroutine to remove repretitive code
    call create_scale_c2l(bounds, c2l_scale_type, &
         scale_c2l(bounds%begc:bounds%endc))
    
    tarr(bounds%begt : bounds%endt) = spval
    sumwt(bounds%begt : bounds%endt) = 0._r8
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .and. col_pp%wttopounit(c) /= 0._r8) then
          l = col_pp%landunit(c)
          if (carr(c) /= spval .and. scale_c2l(c) /= spval .and. scale_l2t(l) /= spval) then
             t = col_pp%wttopounit(c)
             if (sumwt(t) == 0._r8) tarr(t) = 0._r8
             tarr(t) = tarr(t) + carr(c) * scale_c2l(c) * scale_l2t(l) * col_pp%wttopounit(c)
             sumwt(t) = sumwt(t) + col_pp%wttopounit(c)
          end if
       end if
    end do
    found = .false.
    do t = bounds%begt, bounds%endt
       if (sumwt(t) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = t
       else if (sumwt(t) /= 0._r8) then
          tarr(t) = tarr(t)/sumwt(t)
       end if
    end do
    if (found) then
       !#py write(iulog,*)'c2t_1d error: sumwt is greater than 1.0 at t= ',index
       !#py !#py call endrun(decomp_index=index, elmlevel=namet, msg=errMsg(__FILE__, __LINE__))
    end if

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
    real(r8), intent(out) :: tarr( bounds%begt: , 1: ) ! output gridcell array
    character(len=*), intent(in) :: c2l_scale_type     ! scale factor type for averaging
    character(len=*), intent(in) :: l2t_scale_type     ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,t,l,index                       ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_c2l(bounds%begc:bounds%endc)     ! scale factor
    real(r8) :: scale_l2t(bounds%begl:bounds%endl)     ! scale factor
    real(r8) :: sumwt(bounds%begt:bounds%endt)         ! sum of weights
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    call build_scale_l2t(bounds, l2t_scale_type, &
         scale_l2t(bounds%begl:bounds%endl))

    ! Build scale_c2l: TekluTesfa@PNNL created a new subroutine to remove repretitive code
    call create_scale_c2l(bounds, c2l_scale_type, &
         scale_c2l(bounds%begc:bounds%endc))
    
    tarr(bounds%begt : bounds%endt,:) = spval
    do j = 1,num2d
       sumwt(bounds%begt : bounds%endt) = 0._r8
       do c = bounds%begc,bounds%endc 
          if (col_pp%active(c) .and. col_pp%wttopounit(c) /= 0._r8) then
             l = col_pp%landunit(c)
             if (carr(c,j) /= spval .and. scale_c2l(c) /= spval .and. scale_l2t(l) /= spval) then
                t = col_pp%wttopounit(c)
                if (sumwt(t) == 0._r8) tarr(t,j) = 0._r8
                tarr(t,j) = tarr(t,j) + carr(c,j) * scale_c2l(c) * scale_l2t(l) * col_pp%wttopounit(c)
                sumwt(t) = sumwt(t) + col_pp%wttopounit(c)
             end if
          end if
       end do
       found = .false.
       do t = bounds%begt, bounds%endt
          if (sumwt(t) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = t
          else if (sumwt(t) /= 0._r8) then
             tarr(t,j) = tarr(t,j)/sumwt(t)
          end if
       end do
       if (found) then
          !#py write(iulog,*)'c2t_2d error: sumwt is greater than 1.0 at t= ',index
          !#py !#py call endrun(decomp_index=index, elmlevel=namet, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

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
    real(r8), intent(out) :: tarr( bounds%begt: )  ! output gridcell array
    character(len=*), intent(in) :: l2t_scale_type ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: l,t,index                       ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_l2t(bounds%begl:bounds%endl) ! scale factor
    real(r8) :: sumwt(bounds%begt:bounds%endt)     ! sum of weights
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    call build_scale_l2t(bounds, l2t_scale_type, &
         scale_l2t(bounds%begl:bounds%endl))

    tarr(bounds%begt : bounds%endt) = spval
    sumwt(bounds%begt : bounds%endt) = 0._r8
    do l = bounds%begl,bounds%endl
       if (lun_pp%active(l) .and. lun_pp%wttopounit(l) /= 0._r8) then
          if (larr(l) /= spval .and. scale_l2t(l) /= spval) then
             t = lun_pp%wttopounit(l)
             if (sumwt(t) == 0._r8) tarr(t) = 0._r8
             tarr(t) = tarr(t) + larr(l) * scale_l2t(l) * lun_pp%wttopounit(l)
             sumwt(t) = sumwt(t) + lun_pp%wttopounit(l)
          end if
       end if
    end do
    found = .false.
    do t = bounds%begt, bounds%endt
       if (sumwt(t) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = t
       else if (sumwt(t) /= 0._r8) then
          tarr(t) = tarr(t)/sumwt(t)
       end if
    end do
    if (found) then
       !#py write(iulog,*)'l2t_1d error: sumwt is greater than 1.0 at t= ',index
       !#py !#py call endrun(decomp_index=index, elmlevel=namet, msg=errMsg(__FILE__, __LINE__))
    end if

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
    real(r8), intent(out) :: tarr( bounds%begt: , 1: ) ! output gridcell array
    character(len=*), intent(in) :: l2t_scale_type     ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,t,l,index                         ! indices
    integer  :: max_lu_per_gcell                       ! max landunits per gridcell; on the fly
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_l2t(bounds%begl:bounds%endl)     ! scale factor
    real(r8) :: sumwt(bounds%begt:bounds%endt)         ! sum of weights
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    call build_scale_l2t(bounds, l2t_scale_type, &
         scale_l2t(bounds%begl:bounds%endl))

    tarr(bounds%begt : bounds%endt, :) = spval
    do j = 1,num2d
       sumwt(bounds%begt : bounds%endt) = 0._r8
       do l = bounds%begl,bounds%endl
          if (lun_pp%active(l) .and. lun_pp%wttopounit(l) /= 0._r8) then
             if (larr(l,j) /= spval .and. scale_l2t(l) /= spval) then
                t = lun_pp%wttopounit(l)
                if (sumwt(t) == 0._r8) tarr(t,j) = 0._r8
                tarr(t,j) = tarr(t,j) + larr(l,j) * scale_l2t(l) * lun_pp%wttopounit(l)
                sumwt(t) = sumwt(t) + lun_pp%wttopounit(l)
             end if
          end if
       end do
       found = .false.
       do t = bounds%begt,bounds%endt
          if (sumwt(t) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index= t
          else if (sumwt(t) /= 0._r8) then
             tarr(t,j) = tarr(t,j)/sumwt(t)
          end if
       end do
       if (found) then
          !#py write(iulog,*)'l2t_2d error: sumwt is greater than 1.0 at t= ',index,' lev= ',j
          !#py !#py call endrun(decomp_index=index, elmlevel=namet, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine l2t_2d

  !-----------------------------------------------------------------------
  subroutine create_scale_l2g_lookup_gpu(l2g_scale_type, scale_lookup)
    !$acc routine seq
    ! DESCRIPTION:
    ! Create a lookup array, scale_lookup(1..max_lunit), which gives the scale factor for
    ! each landunit type depending on l2g_scale_type
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop, istice, istice_mec, istdlak
    use landunit_varcon, only : isturb_MIN, isturb_MAX, max_lunit
    !
    ! !ARGUMENTS:
    !character(len=*), intent(in)  :: l2g_scale_type           ! scale factor type for averaging
    integer    , intent(in) :: l2g_scale_type  !unity =0, natveg = 3, veg =4, ice=5, nonurb=6, lake=7
    real(r8)   , intent(out) :: scale_lookup(max_lunit)  ! scale factor for each landunit type
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


     ! Initialize scale_lookup to spval for all landunits. Thus, any landunit that keeps
     ! the default value will be excluded from grid cell averages.
     scale_lookup(:) = spval

     if (l2g_scale_type == unity) then
        scale_lookup(:) = 1.0_r8
     else if (l2g_scale_type == natveg) then
        scale_lookup(istsoil) = 1.0_r8
     else if (l2g_scale_type == veg) then
        scale_lookup(istsoil) = 1.0_r8
        scale_lookup(istcrop) = 1.0_r8
     else if (l2g_scale_type == ice) then
        scale_lookup(istice) = 1.0_r8
        scale_lookup(istice_mec) = 1.0_r8
     else if (l2g_scale_type == nonurb) then
        scale_lookup(:) = 1.0_r8
        scale_lookup(isturb_MIN:isturb_MAX) = spval
     else if (l2g_scale_type == lake) then
        scale_lookup(istdlak) = 1.0_r8
     else
        !write(iulog,*)'scale_l2g_lookup_array error: scale type ',l2g_scale_type,' not supported'
        !call endrun(msg=errMsg(__FILE__, __LINE__))
     end if

  end subroutine create_scale_l2g_lookup_gpu

  !-----------------------------------------------------------------------
  subroutine t2g_1d(bounds, tarr, garr, t2g_scale_type)
    !
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from topounits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: tarr( bounds%begt: )  ! input topounit array
    real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
    character(len=*), intent(in) :: t2g_scale_type ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: t,g,index                       ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_t2g(bounds%begt:bounds%endt) ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------
    ! for now, assume that this scale type is always 'unity'
    do t = bounds%begt,bounds%endt
       scale_t2g(t) = 1.0_r8
    end do

    garr(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8
    do t = bounds%begt,bounds%endt
       if (top_pp%wtgcell(t) /= 0._r8) then
          if (tarr(t) /= spval .and. scale_t2g(t) /= spval) then
             g = top_pp%gridcell(t)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + tarr(t) * scale_t2g(t) * top_pp%wtgcell(t)
             sumwt(g) = sumwt(g) + top_pp%wtgcell(t)
          end if
       end if
    end do
    found = .false.
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       !#py write(iulog,*)'t2g_1d error: sumwt is greater than 1.0 at g= ',index
       !#py !#py call endrun(decomp_index=index, elmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
    end if

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
    real(r8), intent(in)  :: tarr( bounds%begt: , 1: ) ! input topounit array
    real(r8), intent(out) :: garr( bounds%begg: , 1: ) ! output gridcell array
    character(len=*), intent(in) :: t2g_scale_type     ! scale factor type for averaging
    !
    ! !LOCAL VARIABLES:
    integer  :: j,g,t,index                         ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_t2g(bounds%begt:bounds%endt)     ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)         ! sum of weights
    !------------------------------------------------------------------------
    ! for now, assume that this scale type is always 'unity'
    do t = bounds%begt,bounds%endt
       scale_t2g(t) = 1.0_r8
    end do

    garr(bounds%begg : bounds%endg, :) = spval
    do j = 1,num2d
       sumwt(bounds%begg : bounds%endg) = 0._r8
       do t = bounds%begt,bounds%endt
          if (top_pp%wtgcell(t) /= 0._r8) then
             if (tarr(t,j) /= spval .and. scale_t2g(t) /= spval) then
                g = top_pp%gridcell(t)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + tarr(t,j) * scale_t2g(t) * top_pp%wtgcell(t)
                sumwt(g) = sumwt(g) + top_pp%wtgcell(t)
             end if
          end if
       end do
       found = .false.
       do g = bounds%begg,bounds%endg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index= g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          !#py write(iulog,*)'t2g_2d error: sumwt is greater than 1.0 at g= ',index,' lev= ',j
          !#py !#py call endrun(decomp_index=index, elmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine t2g_2d

  !-----------------------------------------------------------------------
  subroutine t2g_1d_gpu(bounds, tarr, garr)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from topounits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: tarr( bounds%begt: )  ! input topounit array
    real(r8), intent(out) :: garr( bounds%begg: )  ! output gridcell array
    !
    ! !LOCAL VARIABLES:
    integer  :: t,g,index                       ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_t2g(bounds%begt:bounds%endt) ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------
    ! for now, assume that this scale type is always 'unity'
    do t = bounds%begt,bounds%endt
       scale_t2g(t) = 1.0_r8
    end do

    garr(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8
    do t = bounds%begt,bounds%endt
       if (top_pp%wtgcell(t) /= 0._r8) then
          if (tarr(t) /= spval .and. scale_t2g(t) /= spval) then
             g = top_pp%gridcell(t)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + tarr(t) * scale_t2g(t) * top_pp%wtgcell(t)
             sumwt(g) = sumwt(g) + top_pp%wtgcell(t)
          end if
       end if
    end do
    found = .false.
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
      stop
    end if

  end subroutine t2g_1d_gpu

  !-----------------------------------------------------------------------
  subroutine t2g_2d_gpu(bounds, num2d, tarr, garr)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from topounits to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)  :: num2d                     ! size of second dimension
    real(r8), intent(in)  :: tarr( bounds%begt: , 1: ) ! input topounit array
    real(r8), intent(out) :: garr( bounds%begg: , 1: ) ! output gridcell array
    !
    ! !LOCAL VARIABLES:
    integer  :: j,g,t,index                         ! indices
    logical  :: found                                  ! temporary for error check
    real(r8) :: scale_t2g(bounds%begt:bounds%endt)     ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)         ! sum of weights
    !------------------------------------------------------------------------
    ! for now, assume that this scale type is always 'unity'
    do t = bounds%begt,bounds%endt
       scale_t2g(t) = 1.0_r8
    end do

    garr(bounds%begg : bounds%endg, :) = spval
    do j = 1,num2d
       sumwt(bounds%begg : bounds%endg) = 0._r8
       do t = bounds%begt,bounds%endt
          if (top_pp%wtgcell(t) /= 0._r8) then
             if (tarr(t,j) /= spval .and. scale_t2g(t) /= spval) then
                g = top_pp%gridcell(t)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + tarr(t,j) * scale_t2g(t) * top_pp%wtgcell(t)
                sumwt(g) = sumwt(g) + top_pp%wtgcell(t)
             end if
          end if
       end do
       found = .false.
       do g = bounds%begg,bounds%endg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index= g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          print *, 't2g_2d error: sumwt is greater than 1.0 at g= ',index,' lev= ',j
       end if
    end do

  end subroutine t2g_2d_gpu

  subroutine create_scale_c2l(bounds,c2l_scale_type, scale_c2l)
    ! !DESCRIPTION:
    ! Fill the scale_c2l(bounds%begc:bounds%endc) array with appropriate values for the given c2l_scale_type.
    ! This array can later be used to scale each column in forming higher level averages.
    !
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: c2l_scale_type
    real(r8), intent(inout) :: scale_c2l(bounds%begc:bounds%endc)     ! scale factor
    !-------------- local ----------!
    integer :: c,l
    !-------------------------------!
    if (c2l_scale_type == 'unity') then
      do c = bounds%begc,bounds%endc
          scale_c2l(c) = 1.0_r8
      end do
    else if (c2l_scale_type == 'urbanf') then
      do c = bounds%begc,bounds%endc
          l = col_pp%landunit(c)
          if (lun_pp%urbpoi(l)) then
            if (col_pp%itype(c) == icol_sunwall) then
              scale_c2l(c) = 3.0 * lun_pp%canyon_hwr(l)
            else if (col_pp%itype(c) == icol_shadewall) then
              scale_c2l(c) = 3.0 * lun_pp%canyon_hwr(l)
            else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
              scale_c2l(c) = 3.0_r8
            else if (col_pp%itype(c) == icol_roof) then
              scale_c2l(c) = 1.0_r8
            end if
          else
            scale_c2l(c) = 1.0_r8
          end if
        end do
      else if (c2l_scale_type == 'urbans') then
        do c = bounds%begc,bounds%endc
          l = col_pp%landunit(c)
          if (lun_pp%urbpoi(l)) then
            if (col_pp%itype(c) == icol_sunwall) then
              scale_c2l(c) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
            else if (col_pp%itype(c) == icol_shadewall) then
              scale_c2l(c) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
            else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
              scale_c2l(c) = 3.0 / (2.*lun_pp%canyon_hwr(l) + 1.)
            else if (col_pp%itype(c) == icol_roof) then
              scale_c2l(c) = 1.0_r8
            end if
          else
            scale_c2l(c) = 1.0_r8
          end if
        end do
      else
        !#py write(iulog,*) 'error: scale type ',c2l_scale_type,' not supported'
        !#py !#py call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end subroutine create_scale_c2l

    subroutine create_scale_c2l_gpu(bounds,c2l_scale_type, scale_c2l)
      ! !DESCRIPTION:
      ! Fill the scale_c2l(bounds%begc:bounds%endc) array with appropriate values for the given c2l_scale_type.
      ! This array can later be used to scale each column in forming higher level averages.
      !
      !$acc routine seq
      type(bounds_type), intent(in) :: bounds
      integer , intent(in) :: c2l_scale_type
      real(r8), intent(inout) :: scale_c2l(bounds%begc:bounds%endc)     ! scale factor
      !-------------- local ----------!
      integer :: c,l
      !-----
      if (c2l_scale_type == unity) then
        do c = bounds%begc,bounds%endc
            scale_c2l(c) = 1.0_r8
        end do
      else if (c2l_scale_type == urbanf) then
        do c = bounds%begc,bounds%endc
            l = col_pp%landunit(c)
            if (lun_pp%urbpoi(l)) then
              if (col_pp%itype(c) == icol_sunwall) then
                scale_c2l(c) = 3.0 * lun_pp%canyon_hwr(l)
              else if (col_pp%itype(c) == icol_shadewall) then
                scale_c2l(c) = 3.0 * lun_pp%canyon_hwr(l)
              else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0_r8
              else if (col_pp%itype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
              end if
            else
              scale_c2l(c) = 1.0_r8
            end if
          end do
        else if (c2l_scale_type == urbans) then
          do c = bounds%begc,bounds%endc
            l = col_pp%landunit(c)
            if (lun_pp%urbpoi(l)) then
              if (col_pp%itype(c) == icol_sunwall) then
                scale_c2l(c) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
              else if (col_pp%itype(c) == icol_shadewall) then
                scale_c2l(c) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
              else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0 / (2.*lun_pp%canyon_hwr(l) + 1.)
              else if (col_pp%itype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
              end if
            else
              scale_c2l(c) = 1.0_r8
            end if
          end do
        else
          print *,'error: scale type ',c2l_scale_type,' not supported'
        end if

      end subroutine create_scale_c2l_gpu

       subroutine initialize_scale_l2g_lookup()
         ! DESCRIPTION: 
         ! create module copy of scale_l2g_lookup to avoid creating it for
         ! every subgridAve call. 
         !
         use landunit_varcon, only : istsoil, istcrop, istice, istice_mec, istdlak
         use landunit_varcon, only : isturb_MIN, isturb_MAX, max_lunit
         
         implicit none 

         integer, parameter :: num_scale_types = 6 ! unity, natveg, veg, ice, nonurb, lake
         integer :: itype, lunit 

         allocate(main_scale_l2g_lookup(num_scale_types, max_lunit))

         !itype = unity 
         !$acc parallel loop independent gang vector async(1) 
         do lunit = 1,max_lunit
            main_scale_l2g_lookup(unity,lunit) = 1.0_r8
         end do 

         !itype = natveg
         !$acc parallel loop independent gang vector async(2) 
         do lunit = 1,max_lunit
            if(lunit == istsoil) then 
               main_scale_l2g_lookup(natveg,lunit) = 1.0_r8
            else 
               main_scale_l2g_lookup(natveg,lunit) = spval
            end if 
         end do 

         !itype = veg 
         !$acc parallel loop independent gang vector async(3) 
         do lunit = 1,max_lunit
            if(lunit == istsoil .or. lunit == istcrop ) then 
               main_scale_l2g_lookup(veg,lunit) = 1.0_r8
            else 
               main_scale_l2g_lookup(veg,lunit) = spval
            end if 
         end do 
         
         !itype = ice 
         !$acc parallel loop independent gang vector async(4) 
         do lunit = 1,max_lunit
            if(lunit == istice .or. lunit == istice_mec ) then 
               main_scale_l2g_lookup(ice,lunit) = 1.0_r8
            else 
               main_scale_l2g_lookup(ice,lunit) = spval
            end if 
         end do 

         !itype = nonurb
         !$acc parallel loop independent gang vector async(5) 
         do lunit = 1,max_lunit
            main_scale_l2g_lookup(ice,lunit) = 1.0_r8
            if(lunit >= isturb_MIN .and. lunit <= isturb_MAX) then  
               main_scale_l2g_lookup(ice,lunit) = spval
            end if 
         end do 

         !itype = lake
         !$acc parallel loop independent gang vector async(6) 
         do lunit = 1,max_lunit
            if(lunit == istdlak) then 
               main_scale_l2g_lookup(ice,lunit) = 1.0_r8
            else 
               main_scale_l2g_lookup(ice,lunit) = spval
            end if 
         end do 
         !$acc wait 


       end subroutine initialize_scale_l2g_lookup

       subroutine initialize_scale_c2l(bounds)
         ! Description:
         ! intialize and create scale_c2l to avoid having to calculate it 
         ! every time a c2g or similar routine is called. 
         ! May have to be called again if lun_pp%canyon_hwr is modified.
         !
         ! Note: Could save memory by just allocating over landunits?
         ! Input variables: 
         !
         type(bounds_type), intent(in) :: bounds ! processor level bounds!
         ! Local Variables
         integer :: c, l, begc,endc 
         integer, parameter :: num_scale_types = 3 ! unity, urbanf, urbans 

         begc = bounds%begc
         endc = bounds%endc 
         ! allocate memory :
         allocate(main_scale_c2l(begc:endc,num_scale_types))
         
         ! c2l_scale_type == unity
         !$acc parallel loop independent gang vector default(present) async(1)
         do c = begc,endc
            main_scale_c2l(c,unity) = 1.0_r8
         end do
         
         !c2l_scale_type == urbanf
         !$acc parallel loop independent gang vector default(present) async(2)
         do c = begc,endc
            l = col_pp%landunit(c)
            if (lun_pp%urbpoi(l)) then
               if (col_pp%itype(c) == icol_sunwall) then
                  main_scale_c2l(c,urbanf) = 3.0 * lun_pp%canyon_hwr(l)
               else if (col_pp%itype(c) == icol_shadewall) then
                  main_scale_c2l(c,urbanf) = 3.0 * lun_pp%canyon_hwr(l)
               else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
                  main_scale_c2l(c,urbanf) = 3.0_r8
               else if (col_pp%itype(c) == icol_roof) then
                  main_scale_c2l(c,urbanf) = 1.0_r8
               end if
            else
               main_scale_c2l(c,urbanf) = 1.0_r8
            end if
         end do
         
           ! c2l_scale_type == urbans 
         !$acc parallel loop independent gang vector default(present) async(3)
         do c = begc, endc
           l = col_pp%landunit(c)
           if (lun_pp%urbpoi(l)) then
             if (col_pp%itype(c) == icol_sunwall) then
               main_scale_c2l(c,urbans) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
             else if (col_pp%itype(c) == icol_shadewall) then
               main_scale_c2l(c,urbans) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
             else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
               main_scale_c2l(c,urbans) = 3.0 / (2.*lun_pp%canyon_hwr(l) + 1.)
             else if (col_pp%itype(c) == icol_roof) then
               main_scale_c2l(c,urbans) = 1.0_r8
             end if
           else
            main_scale_c2l(c,urbans) = 1.0_r8
           end if
         end do

         !$acc wait 

       end subroutine initialize_scale_c2l
       
end module subgridAveMod
