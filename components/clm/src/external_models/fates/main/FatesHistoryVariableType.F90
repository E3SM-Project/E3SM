module FatesHistoryVariableType

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals, only : fates_log
  use FatesIOVariableKindMod, only : fates_io_variable_kind_type

  implicit none

  ! This type is instanteated in the HLM-FATES interface (clmfates_interfaceMod.F90)

  type fates_history_variable_type
     character(len=32)    :: vname
     character(len=24)    :: units
     character(len=128)   :: long
     character(len=24)    :: use_default ! States whether a variable should be turned
                                         ! on the output files by default (active/inactive)
                                         ! It is a good idea to set inactive for very large
                                         ! or infrequently used output datasets
     character(len=24)    :: vtype
     character(len=1)     :: avgflag
     integer              :: upfreq  ! Update frequency (this is for checks and flushing)
                                     ! 1 = dynamics "dyn" (daily)
                                     ! 2 = production "prod" (prob model tstep)
     real(r8)             :: flushval
     integer :: dim_kinds_index
     ! Pointers (only one of these is allocated per variable)
     real(r8), pointer     :: r81d(:)
     real(r8), pointer     :: r82d(:,:)
     real(r8), pointer     :: r83d(:,:,:)
     integer,  pointer     :: int1d(:)
     integer,  pointer     :: int2d(:,:)
     integer,  pointer     :: int3d(:,:,:)
   contains
     procedure, public :: Init
     procedure, public :: Flush
     procedure, private :: GetBounds
  end type fates_history_variable_type

contains

  subroutine Init(this, vname, units, long, use_default, &
       vtype, avgflag, flushval, upfreq, num_dim_kinds, dim_kinds, dim_bounds)

    use FatesIODimensionsMod, only : fates_io_dimension_type

    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
    use FatesIOVariableKindMod, only : iotype_index

    implicit none

    class(fates_history_variable_type), intent(inout) :: this
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: long
    character(len=*), intent(in) :: use_default
    character(len=*), intent(in) :: vtype
    character(len=*), intent(in) :: avgflag
    real(r8), intent(in) :: flushval ! If the type is an int we will round with nint
    integer, intent(in) :: upfreq
    integer, intent(in) :: num_dim_kinds
    type(fates_io_dimension_type), intent(in) :: dim_bounds(:)
    type(fates_io_variable_kind_type), intent(inout) :: dim_kinds(:)

    integer :: dk_index
    integer :: lb1, ub1, lb2, ub2
    
    this%vname = vname
    this%units = units
    this%long  = long
    this%use_default = use_default
    this%vtype = vtype
    this%avgflag = avgflag
    this%flushval = flushval
    this%upfreq = upfreq

    nullify(this%r81d)
    nullify(this%r82d)
    nullify(this%r83d)
    nullify(this%int1d)
    nullify(this%int2d)
    nullify(this%int3d)

    dk_index = iotype_index(trim(vtype), num_dim_kinds, dim_kinds)
    this%dim_kinds_index = dk_index
    call dim_kinds(dk_index)%set_active()
                
    call this%GetBounds(0, dim_bounds, dim_kinds, lb1, ub1, lb2, ub2)

    ! NOTE(rgk, 2016-09) currently, all array spaces are flushed each
    ! time the update is called. The flush here on the initialization
    ! may be redundant, but will prevent issues in the future if we
    ! have host models where not all threads are updating the HHistory
    ! array spaces.

    select case(trim(vtype))
    case(patch_r8)
       allocate(this%r81d(lb1:ub1))
       this%r81d(:) = flushval

    case(site_r8)
       allocate(this%r81d(lb1:ub1))
       this%r81d(:) = flushval

    case(patch_ground_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(patch_size_pft_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_ground_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_size_pft_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_size_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_pft_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_age_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_fuel_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_cwdsc_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_can_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_cnlf_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_cnlfpft_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case(site_scag_r8)
       allocate(this%r82d(lb1:ub1, lb2:ub2))
       this%r82d(:,:) = flushval

    case default
       write(fates_log(),*) 'Incompatible vtype passed to set_history_var'
       write(fates_log(),*) 'vtype = ',trim(vtype),' ?'
       stop
       ! end_run
    end select
    
  end subroutine Init
  
  ! =====================================================================================
        
  subroutine GetBounds(this, thread, dim_bounds, dim_kinds, lb1, ub1, lb2, ub2)

    use FatesIODimensionsMod, only : fates_io_dimension_type

    implicit none

     class(fates_history_variable_type), intent(inout) :: this
     integer, intent(in)  :: thread
     type(fates_io_dimension_type), intent(in) :: dim_bounds(:)
     type(fates_io_variable_kind_type), intent(in) :: dim_kinds(:)
     integer, intent(out) :: lb1
     integer, intent(out) :: ub1
     integer, intent(out) :: lb2
     integer, intent(out) :: ub2

     ! local
     integer :: ndims
     integer :: d_index

     lb1 = 0
     ub1 = 0
     lb2 = 0
     ub2 = 0

     ndims = dim_kinds(this%dim_kinds_index)%ndims

     ! The thread = 0 case is the boundaries for the whole proc/node
     if (thread==0) then
        d_index = dim_kinds(this%dim_kinds_index)%dim1_index
        lb1 = dim_bounds(d_index)%lower_bound
        ub1 = dim_bounds(d_index)%upper_bound
        if(ndims>1)then
           d_index = dim_kinds(this%dim_kinds_index)%dim2_index
           lb2 = dim_bounds(d_index)%lower_bound
           ub2 = dim_bounds(d_index)%upper_bound
        end if
     else
        
        d_index = dim_kinds(this%dim_kinds_index)%dim1_index
        lb1 = dim_bounds(d_index)%clump_lower_bound(thread)
        ub1 = dim_bounds(d_index)%clump_upper_bound(thread)
        
        if(ndims>1)then
           d_index = dim_kinds(this%dim_kinds_index)%dim2_index
           lb2 = dim_bounds(d_index)%clump_lower_bound(thread)
           ub2 = dim_bounds(d_index)%clump_upper_bound(thread)
        end if
        
     end if
     
   end subroutine GetBounds

   subroutine Flush(this, thread, dim_bounds, dim_kinds)

    use FatesIODimensionsMod, only : fates_io_dimension_type
    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8, patch_int
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8

    implicit none

    class(fates_history_variable_type), intent(inout) :: this
    integer, intent(in) :: thread
    type(fates_io_dimension_type), intent(in) :: dim_bounds(:)
    type(fates_io_variable_kind_type), intent(in) :: dim_kinds(:)

    integer :: lb1, ub1, lb2, ub2
    
    call this%GetBounds(thread, dim_bounds, dim_kinds, lb1, ub1, lb2, ub2)

    select case(trim(dim_kinds(this%dim_kinds_index)%name))
    case(patch_r8) 
       this%r81d(lb1:ub1) = this%flushval
    case(site_r8) 
       this%r81d(lb1:ub1) = this%flushval
    case(patch_ground_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(patch_size_pft_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_ground_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_size_pft_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_size_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_pft_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_age_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_fuel_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_cwdsc_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_can_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_cnlf_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_cnlfpft_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(site_scag_r8) 
       this%r82d(lb1:ub1, lb2:ub2) = this%flushval
    case(patch_int)
       this%int1d(lb1:ub1) = nint(this%flushval)
    case default
       write(fates_log(),*) 'fates history variable type undefined while flushing history variables'
       stop
       !end_run
    end select
    
 end subroutine Flush

end module FatesHistoryVariableType
