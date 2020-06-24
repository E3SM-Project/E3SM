module histGPUMod
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use decompMod     , only : bounds_type
  use histFileMod   , only : clmptr_rs, clmptr_ra
  use clm_varcon    , only : spval, ispval, dzsoi_decomp
  use clm_varcon    , only : grlnd, nameg, namet, namel, namec, namep
  use GridcellType  , only : grc_pp
  use LandunitType  , only : lun_pp
  use ColumnType    , only : col_pp
  use VegetationType, only : veg_pp

  implicit None

  integer , private, parameter :: no_snow_MIN = 1                 ! minimum valid value for this flag
  integer , public , parameter :: no_snow_normal = 1              ! normal treatment, which should be used for most fields (use spval when snow layer not present)
  integer , public , parameter :: no_snow_zero = 2                ! average in a 0 value for times when the snow layer isn't present
  integer , private, parameter :: no_snow_MAX = 2                 ! maximum valid value for this flag
  integer , private, parameter :: no_snow_unset = no_snow_MIN - 1 ! flag specifying that field is NOT a multi-layer snow field
  integer, parameter :: unity = 0, urbanf = 1, urbans = 2
  integer, parameter :: natveg = 3, veg =4, ice=5, nonurb=6, lake=7
  PUBLIC

  type, public :: history_entry_gpu
     character(len=64), pointer :: name       => null() ! field name
     character(len=16), pointer :: type1d     => null() ! pointer to first dimension type from data type (nameg, etc)
     character(len=16), pointer :: type1d_out => null() ! hbuf first dimension type from data type (nameg, etc)
     character(len=16), pointer :: type2d     => null() ! hbuf second dimension type ["levgrnd","levlak","numrad","ltype","natpft","cft","glc_nec","elevclas","subname(n)","month"]
     integer, pointer :: beg1d   => null()   ! on-node 1d clm pointer start index
     integer, pointer :: end1d   => null()   ! on-node 1d clm pointer end index
     integer, pointer :: num1d   => null()   ! size of clm pointer first dimension (all nodes)
     integer, pointer :: numdims => null()   ! the actual number of dimensions, this allows
                                             ! for 2D arrays, where the second dimension is allowed
                                             ! to be 1
     integer,pointer :: beg1d_out =>null()   ! on-node 1d hbuf pointer start index
     integer,pointer :: end1d_out =>null()   ! on-node 1d hbuf pointer end index
     integer,pointer :: num1d_out =>null()   ! size of hbuf first dimension (all nodes)
     integer,pointer :: num2d     =>null()   ! size of hbuf second dimension (e.g. number of vertical levels)
     integer,pointer :: hpindex   =>null()   ! history pointer index
     integer, pointer :: p2c_scale_type => null() ! scale factor when averaging pft to column
     integer, pointer :: c2l_scale_type => null() ! scale factor when averaging column to landunit
     integer, pointer :: l2g_scale_type => null() ! scale factor when averaging landunit to gridcell
     character(len=8), pointer :: t2g_scale_type => null() ! scale factor when averaging topounit to gridcell
     integer, pointer :: no_snow_behavior => null()        ! for multi-layer snow fields, flag saying how to treat times when a given snow layer is absent
     !
     character(len=1), pointer  :: avgflag => null()   ! time averaging flag
     real(r8), pointer :: hbuf(:,:)   => null()         ! history buffer (dimensions: dim1d x num2d)
     integer , pointer :: nacs(:,:)   => null()         ! accumulation counter (dimensions: dim1d x num2d)
  end type history_entry_gpu

  type, public :: history_tape_gpu
     integer , pointer  :: nflds      => null() ! number of active fields on tape
     integer , pointer  :: ntimes     => null() ! current number of time samples on tape
     integer , pointer  :: mfilt      => null() ! maximum number of time samples per tape
     integer , pointer  :: nhtfrq     => null() ! number of time samples per tape
     logical , pointer  :: is_endhist => null() ! true => current time step is end of history interval
     real(r8), pointer  :: begtime    => null() ! time at beginning of history averaging interval
     type (history_entry_gpu), allocatable :: hlist(:)   ! array of active history tape entries
  end type history_tape_gpu

  type (history_tape_gpu), public, allocatable :: tape_gpu(:)    ! array history tapes
  !$acc declare create(tape_gpu)

contains

  subroutine htape_gpu_init()
    use histFileMod, only : tape, ntapes

    implicit none
    integer :: size1,size2,t,f
    ! ---- subroutine to store tape data into type for gpu
    print *, "allocating gpu ",ntapes,"tapes"
    allocate(tape_gpu(ntapes))
    do t = 1, ntapes
      print *, "allocating hlist with ",tape(t)%nflds,"fields"
      allocate(tape_gpu(t)%hlist(tape(t)%nflds))
      allocate(tape_gpu(t)%nflds )
      allocate(tape_gpu(t)%ntimes)
      allocate(tape_gpu(t)%mfilt )
      allocate(tape_gpu(t)%nhtfrq)
      tape_gpu(t)%nflds  = tape(t)%nflds
      tape_gpu(t)%ntimes = tape(t)%ntimes
      tape_gpu(t)%mfilt  = tape(t)%mfilt
      tape_gpu(t)%nhtfrq = tape(t)%nhtfrq

      do f = 1, tape_gpu(t)%nflds
        size1 = size(tape(t)%hlist(f)%hbuf,1);
        size2 = size(tape(t)%hlist(f)%hbuf,2)
        allocate(tape_gpu(t)%hlist(f)%hbuf(size1,size2))
        size1 = size(tape(t)%hlist(f)%nacs,1);
        size2 = size(tape(t)%hlist(f)%nacs,2)
        allocate(tape_gpu(t)%hlist(f)%nacs(size1,size2))
        tape_gpu(t)%hlist(f)%hbuf(:,:) = tape(t)%hlist(f)%hbuf(:,:)
        tape_gpu(t)%hlist(f)%nacs(:,:) = tape(t)%hlist(f)%nacs(:,:)
        !!
        allocate(tape_gpu(t)%hlist(f)%avgflag)
        tape_gpu(t)%hlist(f)%avgflag = tape(t)%hlist(f)%avgflag
        allocate(tape_gpu(t)%hlist(f)%name            )
        allocate(tape_gpu(t)%hlist(f)%type1d          )
        allocate(tape_gpu(t)%hlist(f)%type1d_out      )
        allocate(tape_gpu(t)%hlist(f)%type2d          )
        allocate(tape_gpu(t)%hlist(f)%beg1d           )
        allocate(tape_gpu(t)%hlist(f)%end1d           )
        allocate(tape_gpu(t)%hlist(f)%num1d           )
        allocate(tape_gpu(t)%hlist(f)%numdims         )
        allocate(tape_gpu(t)%hlist(f)%beg1d_out       )
        allocate(tape_gpu(t)%hlist(f)%end1d_out       )
        allocate(tape_gpu(t)%hlist(f)%num1d_out       )
        allocate(tape_gpu(t)%hlist(f)%num2d           )
        allocate(tape_gpu(t)%hlist(f)%hpindex         )
        allocate(tape_gpu(t)%hlist(f)%p2c_scale_type  )
        allocate(tape_gpu(t)%hlist(f)%c2l_scale_type  )
        allocate(tape_gpu(t)%hlist(f)%l2g_scale_type  )
        allocate(tape_gpu(t)%hlist(f)%t2g_scale_type  )
        allocate(tape_gpu(t)%hlist(f)%no_snow_behavior)
        tape_gpu(t)%hlist(f)%name             = tape(t)%hlist(f)%field%name
        tape_gpu(t)%hlist(f)%type1d           = tape(t)%hlist(f)%field%type1d
        tape_gpu(t)%hlist(f)%type1d_out       = tape(t)%hlist(f)%field%type1d_out
        tape_gpu(t)%hlist(f)%type2d           = tape(t)%hlist(f)%field%type2d
        tape_gpu(t)%hlist(f)%beg1d            = tape(t)%hlist(f)%field%beg1d
        tape_gpu(t)%hlist(f)%end1d            = tape(t)%hlist(f)%field%end1d
        tape_gpu(t)%hlist(f)%num1d            = tape(t)%hlist(f)%field%num1d
        tape_gpu(t)%hlist(f)%numdims          = tape(t)%hlist(f)%field%numdims
        tape_gpu(t)%hlist(f)%beg1d_out        = tape(t)%hlist(f)%field%beg1d_out
        tape_gpu(t)%hlist(f)%end1d_out        = tape(t)%hlist(f)%field%end1d_out
        tape_gpu(t)%hlist(f)%num1d_out        = tape(t)%hlist(f)%field%num1d_out
        tape_gpu(t)%hlist(f)%num2d            = tape(t)%hlist(f)%field%num2d
        tape_gpu(t)%hlist(f)%hpindex          = tape(t)%hlist(f)%field%hpindex
        if(trim(tape(t)%hlist(f)%field%p2c_scale_type) =='unity') then
            tape_gpu(t)%hlist(f)%p2c_scale_type = unity
        elseif(trim(tape(t)%hlist(f)%field%p2c_scale_type) =='urbanf')Then
            tape_gpu(t)%hlist(f)%p2c_scale_type = urbanf
        elseif(trim(tape(t)%hlist(f)%field%p2c_scale_type) =='urbans')Then
            tape_gpu(t)%hlist(f)%p2c_scale_type = urbans
        end if
        if(trim(tape(t)%hlist(f)%field%c2l_scale_type) =='unity') then
            tape_gpu(t)%hlist(f)%c2l_scale_type = unity
        elseif(trim(tape(t)%hlist(f)%field%c2l_scale_type) =='urbanf')Then
            tape_gpu(t)%hlist(f)%c2l_scale_type = urbanf
        elseif(trim(tape(t)%hlist(f)%field%c2l_scale_type) =='urbans')Then
            tape_gpu(t)%hlist(f)%c2l_scale_type = urbans
        end if

        if (trim(tape(t)%hlist(f)%field%l2g_scale_type) == 'unity') then
           tape_gpu(t)%hlist(f)%l2g_scale_type = unity
        else if (trim(tape(t)%hlist(f)%field%l2g_scale_type) == 'natveg') then
           tape_gpu(t)%hlist(f)%l2g_scale_type = natveg
        else if (trim(tape(t)%hlist(f)%field%l2g_scale_type) == 'veg') then
           tape_gpu(t)%hlist(f)%l2g_scale_type = veg
        else if (trim(tape(t)%hlist(f)%field%l2g_scale_type) == 'ice') then
           tape_gpu(t)%hlist(f)%l2g_scale_type = ice
        else if (trim(tape(t)%hlist(f)%field%l2g_scale_type) == 'nonurb') then
           tape_gpu(t)%hlist(f)%l2g_scale_type = nonurb
        else if (trim(tape(t)%hlist(f)%field%l2g_scale_type) == 'lake') then
           tape_gpu(t)%hlist(f)%l2g_scale_type = lake
        else
           print *, "scale_l2g_lookup_array error: scale type  not supported"
        end if
        
        tape_gpu(t)%hlist(f)%t2g_scale_type   = tape(t)%hlist(f)%field%t2g_scale_type
        tape_gpu(t)%hlist(f)%no_snow_behavior = tape(t)%hlist(f)%field%no_snow_behavior
      end do
    enddo

  end subroutine htape_gpu_init

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf_gpu(step,inc, nclumps)
    !
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    use histFileMod, only : ntapes
    use decompMod, only : get_proc_bounds,get_clump_bounds_gpu, bounds_type
    ! !ARGUMENTS:
    integer , intent(in) :: step
    integer , intent(in) :: inc
    integer, value, intent(in) :: nclumps 
    !
    ! !LOCAL VARIABLES:
    type(bounds_type) :: bounds
    integer :: t , hp, nc          ! tape index
    integer :: f                   ! field index
    integer :: numdims             ! number of dimensions
    integer :: num2d               ! size of second dimension (e.g. number of vertical levels)
    !----------------------------------------------------------------------
    !$acc parallel vector_length(32) default(present)
    !$acc loop independent gang  private(nc,t,f,hp,numdims,num2d,bounds)
    do nc = 1, nclumps
            !call get_proc_bounds(bounds) 
            call get_clump_bounds_gpu(nc,bounds)
              !$acc loop seq
             do t = 1,ntapes
                !$acc loop seq
                do f = 1,tape_gpu(t)%nflds
                         numdims = tape_gpu(t)%hlist(f)%numdims
                        
                        if ( numdims == 1) then
                                hp = tape_gpu(t)%hlist(f)%hpindex            
                                call hist_update_hbuf_field_1d_gpu(t, f, hp ,bounds)
                        else
                                hp = tape_gpu(t)%hlist(f)%hpindex
                                num2d = tape_gpu(t)%hlist(f)%num2d
                                call hist_update_hbuf_field_2d_gpu (t, f, hp , bounds, num2d)
                        end if
                end do
          end do
     end do 
    !$acc end parallel
   
    !TODO: change inc to be end of hist interval?   
    if(mod(step,inc) == 0 .and. step .ne. 0) then
        print *, "transfering tape to cpu:"
        call transfer_tape_to_cpu()
    endif
  end subroutine hist_update_hbuf_gpu

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf_field_1d_gpu(t, f, hpindex, bounds)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! This canNOT be called from within a threaded region (see comment below regarding the
    ! call to p2g, and the lack of explicit bounds on its arguments; see also bug 1786)
    !
    ! !USES:
    use subgridAveMod   , only : p2g, c2g, l2g, t2g
    use subgridAveMod   , only : unity, urbanf, urbans, natveg
    use subgridAveMod   , only : natveg, veg,ice,nonurb,lake
    use landunit_varcon , only : istice_mec
    !use decompMod       , only : BOUNDS_LEVEL_PROCdd

    ! !ARGUMENTS:
    integer,value, intent(in) :: t            ! tape index
    integer,value, intent(in) :: f            ! field index
    integer,value, intent(in) :: hpindex
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: k                       ! gridcell, landunit, column or pft index
    logical  :: check_active            ! true => check 'active' flag of each point (this refers to a point being active, NOT a history field being active)
    logical  :: valid                   ! true => history operation is valid
    logical  :: map2gcell               ! true => map clm pointer field to gridcell
    logical , pointer :: active(:)      ! flag saying whether each point is active (used for type1d = landunit/column/pft) (this refers to a point being active, NOT a history field being active)
    real(r8) :: field_gcell(bounds%begg:bounds%endg)  ! gricell level field (used if mapping to gridcell is done)
    integer  :: j
    integer  :: k_offset                    ! offset for mapping sliced subarray pointers when outputting variables in PFT/col vector form
    !-----------------------------------------------------------------------

    associate(&
      avgflag        =>  tape_gpu(t)%hlist(f)%avgflag  ,&
      nacs           =>  tape_gpu(t)%hlist(f)%nacs       ,&
      hbuf           =>  tape_gpu(t)%hlist(f)%hbuf       ,&
      beg1d          =>  tape_gpu(t)%hlist(f)%beg1d     ,&
      end1d          =>  tape_gpu(t)%hlist(f)%end1d      ,&
      type1d         =>  tape_gpu(t)%hlist(f)%type1d     ,&
      type1d_out     =>  tape_gpu(t)%hlist(f)%type1d_out    ,&
      p2c_scale_type =>  tape_gpu(t)%hlist(f)%p2c_scale_type,&
      c2l_scale_type =>  tape_gpu(t)%hlist(f)%c2l_scale_type,&
      l2g_scale_type =>  tape_gpu(t)%hlist(f)%l2g_scale_type,&
      field          =>  clmptr_rs(hpindex)%ptr &
      )
    ! set variables to check weights when allocate all pfts

    map2gcell = .false.
    if (type1d_out == nameg .or. type1d_out == grlnd) then
       if (type1d == namep) then
          ! In this and the following calls, we do NOT explicitly subset field using
          ! bounds (e.g., we do NOT do field(bounds%begp:bounds%endp). This is because,
          ! for some fields, the lower bound has been reset to 1 due to taking a pointer
          ! to an array slice. Thus, this code will NOT work properly if done within a
          ! threaded region! (See also bug 1786)
          !!if(tape_gpu(t)%hlist(f)%name == 'TREFMNAV') then
          !!        print *,"field:",bounds%begp
          !!        print *, field(bounds%begp:bounds%endp)
          !!end if 
          call p2g(bounds, &
               field(bounds%begp:bounds%endp), &
               field_gcell(bounds%begg:bounds%endg), &
               p2c_scale_type, c2l_scale_type, l2g_scale_type)
          
          !!if(tape_gpu(t)%hlist(f)%name == 'TREFMNAV') then
          !!        print *,"field_gcell:",bounds%begg
          !!        print *, field_gcell(bounds%begg:bounds%endg)
          !!end if 
         
          map2gcell = .true.
       else if (type1d == namec) then
          call c2g(bounds, &
               field(bounds%begc:bounds%endc), &
               field_gcell(bounds%begg:bounds%endg), &
               c2l_scale_type, l2g_scale_type)
          
          map2gcell = .true.
       else if (type1d == namel) then
          call l2g(bounds, &
               field(bounds%begl:bounds%endl), &
               field_gcell(bounds%begg:bounds%endg), &
               l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namet) then
          call t2g(bounds, &
               field(bounds%begt:bounds%endt), &
               field_gcell(bounds%begg:bounds%endg))
          map2gcell = .true.
       end if
    end if

    if (map2gcell) then  ! Map to gridcell
       ! note that in this case beg1d = begg and end1d=endg
       select case (avgflag)
       case ('I') ! Instantaneous
          do k = bounds%begg,bounds%endg
             if (field_gcell(k) /= spval) then
                hbuf(k,1) = field_gcell(k)
             else
                hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('A') ! Time average
          do k = bounds%begg,bounds%endg
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                hbuf(k,1) = hbuf(k,1) + field_gcell(k)
                nacs(k,1) = nacs(k,1) + 1
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
          end do
       case ('X') ! Maximum over time
          do k = bounds%begg,bounds%endg
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = -1.e50_r8
                hbuf(k,1) = max( hbuf(k,1), field_gcell(k) )
             else
                hbuf(k,1) = spval
             endif
             nacs(k,1) = 1
          end do
       case ('M') ! Minimum over time
          do k = bounds%begg,bounds%endg
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = +1.e50_r8
                hbuf(k,1) = min( hbuf(k,1), field_gcell(k) )
             else
                hbuf(k,1) = spval
             endif
             nacs(k,1) = 1
          end do
       case default
          print *,' ERROR: invalid time averaging flag ', avgflag
       end select

    else  ! Do not map to gridcell

       ! For data defined on the pft, col, and landunit we need to check if a point is active
       ! to determine whether that point should be assigned spval
       if (type1d == namep) then
          check_active = .true.
          active => veg_pp%active
       else if (type1d == namec) then
          check_active = .true.
          active => col_pp%active
       else if (type1d == namel) then
          check_active = .true.
          active =>lun_pp%active
       else
          check_active = .false.
       end if

       select case (avgflag)
       case ('I') ! Instantaneous
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   hbuf(k,1) = field(k)
                else
                   hbuf(k,1) = spval
                end if
             else
                hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('A') ! Time average
          ! create mappings for array slice pointers (which go from 1 to size(field) rather than beg1d to end1d)
          if ( end1d .eq. ubound(field,1) ) then
             k_offset = 0
          else
             print *, "koffset != 0"
             k_offset = 1 - beg1d
          endif
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k+k_offset) /= spval) then   ! add k_offset
                   if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                   hbuf(k,1) = hbuf(k,1) + field(k+k_offset)   ! add k_offset
                   nacs(k,1) = nacs(k,1) + 1
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
          end do
       case ('X') ! Maximum over time
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   if (nacs(k,1) == 0) hbuf(k,1) = -1.e50_r8
                   hbuf(k,1) = max( hbuf(k,1), field(k) )
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('M') ! Minimum over time
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   if (nacs(k,1) == 0) hbuf(k,1) = +1.e50_r8
                   hbuf(k,1) = min( hbuf(k,1), field(k) )
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case default
          print *,' ERROR: invalid time averaging flag ', avgflag
          !call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    end if
  end associate
end subroutine hist_update_hbuf_field_1d_gpu

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf_field_2d_gpu (t, f,hpindex, bounds, num2d)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! This canNOT be called from within a threaded region (see comment below regarding the
    ! call to p2g, and the lack of explicit bounds on its arguments; see also bug 1786)
    !
    ! !USES:
    use subgridAveMod   , only : p2g, c2g, l2g, t2g
    use landunit_varcon , only : istice_mec
    !use decompMod       , only : BOUNDS_LEVEL_PROC
    !
    ! !ARGUMENTS:
    integer,value, intent(in) :: t            ! tape index
    integer,value, intent(in) :: f            ! field index
    integer,value, intent(in) :: hpindex
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num2d        ! size of second dimension
    !
    ! !LOCAL VARIABLES:
    integer  :: k                       ! gridcell, landunit, column or pft index
    integer  :: j                       ! level index
    integer  :: beg1d,end1d             ! beginning and ending indices
    logical  :: check_active            ! true => check 'active' flag of each point (this refers to a point being active, NOT a history field being active)
    logical  :: valid                   ! true => history operation is valid
    logical  :: map2gcell               ! true => map clm pointer field to gridcell
    logical           :: field_allocated! whether 'field' was allocated here
    logical , pointer :: active(:)      ! flag saying whether each point is active (used for type1d = landunit/column/pft)
                                        !(this refers to a point being active, NOT a history field being active)
    real(r8) :: field_gcell(bounds%begg:bounds%endg,num2d) ! gricell level field (used if mapping to gridcell is done)
    !-----------------------------------------------------------------------

    associate(&
    avgflag             =>  tape_gpu(t)%hlist(f)%avgflag         ,&
    nacs                =>  tape_gpu(t)%hlist(f)%nacs             ,&
    hbuf                =>  tape_gpu(t)%hlist(f)%hbuf             ,&
    beg1d               =>  tape_gpu(t)%hlist(f)%beg1d           ,&
    end1d               =>  tape_gpu(t)%hlist(f)%end1d           ,&
    type1d              =>  tape_gpu(t)%hlist(f)%type1d          ,&
    type1d_out          =>  tape_gpu(t)%hlist(f)%type1d_out      ,&
    p2c_scale_type      =>  tape_gpu(t)%hlist(f)%p2c_scale_type  ,&
    c2l_scale_type      =>  tape_gpu(t)%hlist(f)%c2l_scale_type  ,&
    l2g_scale_type      =>  tape_gpu(t)%hlist(f)%l2g_scale_type  ,&
    no_snow_behavior    =>  tape_gpu(t)%hlist(f)%no_snow_behavior,&
    field               =>  clmptr_ra(hpindex)%ptr &
    )

    if (no_snow_behavior /= no_snow_unset) then
       ! For multi-layer snow fields, build a special output variable that handles
       ! missing snow layers appropriately

       ! Note, regarding bug 1786: The following allocation is not what we would want if
       ! this routine were operating in a threaded region (or, more generally, within a
       ! loop over nclumps) - in that case we would want to use the bounds information for
       ! this clump. But currently that's not possible because the bounds of some fields
       ! have been reset to 1 - see also bug 1786. Similarly, if we wanted to allow
       ! operation within a loop over clumps, we would need to pass 'bounds' to
       ! hist_set_snow_field_2d rather than relying on beg1d & end1d (which give the proc,
       ! bounds not the clump bounds)

       !allocate(field(lbound(clmptr_ra(hpindex)%ptr, 1) : ubound(clmptr_ra(hpindex)%ptr, 1), num2d))
       field_allocated = .true.
       print *, "calling hist set snow field 2d:"
       stop
       !call hist_set_snow_field_2d(field, clmptr_ra(hpindex)%ptr, no_snow_behavior, type1d, &
       !     beg1d, end1d)
    else

       field_allocated = .false.
    end if

    ! set variables to check weights when allocate all pfts

    map2gcell = .false.
    if (type1d_out == nameg .or. type1d_out == grlnd) then
       if (type1d == namep) then
          ! In this and the following calls, we do NOT explicitly subset field using
                ! (e.g., we do NOT do field(bounds%begp:bounds%endp). This is because,
          ! for some fields, the lower bound has been reset to 1 due to taking a pointer
          ! to an array slice. Thus, this code will NOT work properly if done within a
          ! threaded region! (See also bug 1786)
          call p2g(bounds, num2d, &
               field(bounds%begp:bounds%endp,:), &
               field_gcell(bounds%begg:bounds%endg, :), &
               p2c_scale_type, c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namec) then
          call c2g(bounds, num2d, &
               field(bounds%begc:bounds%endc,:), &
               field_gcell(bounds%begg:bounds%endg, :), &
               c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namel) then
          call l2g(bounds, num2d, &
               field(bounds%begl:bounds%endl,:), &
               field_gcell(bounds%begg:bounds%endg, :), &
               l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namet) then
          call t2g(bounds, num2d, &
               field(bounds%begt:bounds%endt,:), &
               field_gcell(bounds%begg:bounds%endg, :))
          map2gcell = .true.
       end if
    end if

    if (map2gcell) then  ! Map to gridcell

       ! note that in this case beg1d = begg and end1d=endg
       select case (avgflag)
       case ('I') ! Instantaneous
          do j = 1,num2d
             do k = bounds%begg,bounds%endg
                if (field_gcell(k,j) /= spval) then
                   hbuf(k,j) = field_gcell(k,j)
                else
                   hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('A') ! Time average
          do j = 1,num2d
             do k = bounds%begg,bounds%endg
                if (field_gcell(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                   hbuf(k,j) = hbuf(k,j) + field_gcell(k,j)
                   nacs(k,j) = nacs(k,j) + 1
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                endif
             end do
          end do
       case ('X') ! Maximum over time
          do j = 1,num2d
             do k = bounds%begg,bounds%endg
                if (field_gcell(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = -1.e50_r8
                   hbuf(k,j) = max( hbuf(k,j), field_gcell(k,j) )
                else
                   hbuf(k,j) = spval
                endif
                nacs(k,j) = 1
             end do
          end do
       case ('M') ! Minimum over time
          do j = 1,num2d
             do k = bounds%begg,bounds%endg
                if (field_gcell(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                   hbuf(k,j) = min( hbuf(k,j), field_gcell(k,j) )
                else
                   hbuf(k,j) = spval
                endif
                nacs(k,j) = 1
             end do
          end do
       case default
          print *,' ERROR: invalid time averaging flag ', avgflag
       end select

    else  ! Do not map to gridcell

       ! For data defined on the pft, col or landunit, we need to check if a point is active
       ! to determine whether that point should be assigned spval
       if (type1d == namep) then
          check_active = .true.
          active => veg_pp%active
       else if (type1d == namec) then
          check_active = .true.
          active => col_pp%active
       else if (type1d == namel) then
          check_active = .true.
          active =>lun_pp%active
       else
          check_active = .false.
       end if

       ! Note that since field points to an array section the
       ! bounds are field(1:end1d-beg1d+1, num2d) - therefore
       ! need to do the shifting below

       select case (avgflag)
       case ('I') ! Instantaneous
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      hbuf(k,j) = field(k-beg1d+1,j)
                   else
                      hbuf(k,j) = spval
                   end if
                else
                   hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('A') ! Time average
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                      hbuf(k,j) = hbuf(k,j) + field(k-beg1d+1,j)
                      nacs(k,j) = nacs(k,j) + 1
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
             end do
          end do
       case ('X') ! Maximum over time
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = -1.e50_r8
                      hbuf(k,j) = max( hbuf(k,j), field(k-beg1d+1,j) )
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('M') ! Minimum over time
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                      hbuf(k,j) = min( hbuf(k,j), field(k-beg1d+1,j))
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case default
          print *,' ERROR: invalid time averaging flag ', avgflag
          !call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    end if

    if (field_allocated) then
       deallocate(field)
    end if
  end associate
  end subroutine hist_update_hbuf_field_2d_gpu

  !-----------------------------------------------------------------------
  subroutine hist_set_snow_field_2d (field_out, field_in, no_snow_behavior, type1d, beg1d, end1d)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Set values in history field dimensioned by levsno.
    !
    ! This routine handles what to do when a given snow layer doesn't exist for a given
    ! point, based on the no_snow_behavior argument. Options are:
    !
    ! - no_snow_normal: This is the normal behavior, which applies to most snow fields:
    !   Use spval (missing value flag). This means that temporal averages will just
    !   consider times when a particular snow layer actually existed
    !
    ! - no_snow_zero: Average in a 0 value for times when the snow layer isn't present
    !
    ! Input and output fields can be defined at the pft or column level
    !
    ! !ARGUMENTS:
    integer         , intent(in)  :: beg1d                    ! beginning spatial index
    integer         , intent(in)  :: end1d                    ! ending spatial index
    real(r8)        , intent(out) :: field_out( beg1d: , 1: ) ! output field [point, lev]
    real(r8)        , intent(in)  :: field_in ( beg1d: , 1: ) ! input field [point, lev]
    integer         , intent(in)  :: no_snow_behavior         ! behavior to use when a snow layer is absent
    character(len=16), intent(in)  :: type1d                   ! 1d clm pointer type ("column" or "pft")
    !
    ! !LOCAL VARIABLES:
    integer :: num_levels             ! total number of possible snow layers
    integer :: point
    integer :: level
    integer :: num_snow_layers        ! number of snow layers that exist at a point
    integer :: num_nonexistent_layers
    integer :: c                      ! column index
    real(r8):: no_snow_val            ! value to use when a snow layer is missing
    !-----------------------------------------------------------------------

    associate(&
    snl            => col_pp%snl  &   ! Input: [integer (:)] number of snow layers (negative)
    )

    num_levels = ubound(field_in, 2)

    ! Determine no_snow_val
    select case (no_snow_behavior)
    case (no_snow_normal)
       no_snow_val = spval
    case (no_snow_zero)
       no_snow_val = 0._r8
    case default
       print *, ' ERROR: unrecognized no_snow_behavior: ', no_snow_behavior
    end select

    do point = beg1d, end1d

       ! Get number of snow layers at this point

       if (type1d == namec) then
          c = point
       else if (type1d == namep) then
          c = veg_pp%column(point)
       else
          print *, ' ERROR: Only implemented for pft and col-level fields'
       end if

       num_snow_layers = abs(snl(c))
       num_nonexistent_layers = num_levels - num_snow_layers

       ! Fill output field appropriately for each layer
       ! When only a subset of snow layers exist, it is the LAST num_snow_layers that exist

       do level = 1, num_nonexistent_layers
          field_out(point, level) = no_snow_val
       end do
       do level = (num_nonexistent_layers + 1), num_levels
          field_out(point, level) = field_in(point, level)
       end do

    end do

    end associate

  end subroutine hist_set_snow_field_2d


  subroutine transfer_tape_to_cpu()
    use histFileMod, only : tape, ntapes

    implicit none
    integer :: size1,size2,t,f

    print *, "updating tape_gpu on cpu:"
    !$acc update self(tape_gpu)

    !loop is done on cpu --- could accelerate using openACC cpu threading?
    print *, "update tape on cpu"
    do t = 1, ntapes
      do f = 1, tape_gpu(t)%nflds
        tape(t)%hlist(f)%hbuf(:,:) = tape_gpu(t)%hlist(f)%hbuf(:,:)
        tape(t)%hlist(f)%nacs(:,:) = tape_gpu(t)%hlist(f)%nacs(:,:)
      end do
    end do


  end subroutine
  
  subroutine set_gpu_tape()
         use histFileMod, only : tape, ntapes
          implicit none
          integer  :: t 
          integer :: f

          do t = 1 , ntapes
              if(tape(t)%is_endhist) then 
                print *,  "adjusting gpu tape after normalization/zeroing",t
                do f = 1, tape_gpu(t)%nflds 
                        tape_gpu(t)%hlist(f)%hbuf(:,:) = 0d0 
                        tape_gpu(t)%hlist(f)%nacs(:,:) = 0
                 end do
                !$acc update device(tape_gpu(t))
                end if  
          end do 

  end subroutine set_gpu_tape


end module histGPUMod
