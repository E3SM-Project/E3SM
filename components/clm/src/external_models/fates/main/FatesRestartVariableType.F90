module FatesRestartVariableMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals, only : fates_log
  use FatesIOVariableKindMod, only : fates_io_variable_kind_type

  implicit none

  ! This type is instanteated in the HLM-FATES interface (clmfates_interfaceMod.F90)
  
  type fates_restart_variable_type
     character(len=32)    :: vname
     character(len=24)    :: units
     character(len=128)   :: long
     character(len=24)    :: vtype
     real(r8)             :: flushval  ! DONT THINK THIS IS NEEDED IN RESTARTS
                                       ! RESTARTS HAVE A MAPPING TABLE AND
                                       ! THERE IS NO AVERAGING AND NO NEED TO
                                       ! INDICATE NON-INCLUDED ARRAY SPACES
                                       ! KEEPING FOR NOW (RGK-11-2016)
     integer :: dim_kinds_index
     ! Pointers (only one of these is allocated per variable)
     real(r8), pointer     :: r81d(:)
     integer,  pointer     :: int1d(:)
   contains
     procedure, public :: Init
     procedure, public :: Flush
     procedure, private :: GetBounds
  end type fates_restart_variable_type

contains

  subroutine Init(this, vname, units, long, vtype, flushval, num_dim_kinds, dim_kinds, dim_bounds)

    use FatesIODimensionsMod, only : fates_io_dimension_type
    use FatesIOVariableKindMod, only : patch_r8, site_r8, cohort_r8
    use FatesIOVariableKindMod, only : patch_int, site_int, cohort_int
    use FatesIOVariableKindMod, only : iotype_index

    implicit none
    
    class(fates_restart_variable_type), intent(inout) :: this
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: long
    character(len=*), intent(in) :: vtype
    real(r8), intent(in) :: flushval
    integer, intent(in) :: num_dim_kinds
    type(fates_io_dimension_type), intent(in) :: dim_bounds(:)
    type(fates_io_variable_kind_type), intent(inout) :: dim_kinds(:)

    integer :: dk_index
    integer :: lb1, ub1, lb2, ub2
    
    this%vname = vname
    this%units = units
    this%long  = long
    this%vtype = vtype
    this%flushval = flushval

    nullify(this%r81d)
    nullify(this%int1d)

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

    case(cohort_r8)
       allocate(this%r81d(lb1:ub1))
       this%r81d(:) = flushval

    case(patch_r8)
       allocate(this%r81d(lb1:ub1))
       this%r81d(:) = flushval

    case(site_r8)
       allocate(this%r81d(lb1:ub1))
       this%r81d(:) = flushval

    case(cohort_int)
       allocate(this%int1d(lb1:ub1))
       this%int1d(:) = idnint(flushval)

    case(patch_int)
       allocate(this%int1d(lb1:ub1))
       this%int1d(:) = idnint(flushval)

    case(site_int)
       allocate(this%int1d(lb1:ub1))
       this%int1d(:) = idnint(flushval)

    case default
       write(fates_log(),*) 'Incompatible vtype passed to set_restart_var'
       write(fates_log(),*) 'vtype = ',trim(vtype),' ?'
       stop
       ! end_run
    end select
    
  end subroutine Init
  
  ! =====================================================================================

  subroutine GetBounds(this, thread, dim_bounds, dim_kinds, lb1, ub1, lb2, ub2)

    use FatesIODimensionsMod, only : fates_io_dimension_type

    implicit none
    
     class(fates_restart_variable_type), intent(inout) :: this
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

   ! ====================================================================================

  subroutine flush(this, thread, dim_bounds, dim_kinds)

    use FatesIODimensionsMod, only : fates_io_dimension_type
    use FatesIOVariableKindMod, only : patch_r8, site_r8, cohort_r8
    use FatesIOVariableKindMod, only : patch_int, site_int, cohort_int

    implicit none

    class(fates_restart_variable_type), intent(inout) :: this
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
    case(cohort_r8)
       this%r81d(lb1:ub1) = this%flushval
    case(patch_int)
       this%int1d(lb1:ub1) = nint(this%flushval)
    case(site_int)
       this%int1d(lb1:ub1) = nint(this%flushval)
    case(cohort_int)
       this%int1d(lb1:ub1) = nint(this%flushval)
       
    case default
       write(fates_log(),*) 'fates history variable type undefined while flushing history variables'
       stop
       !end_run
    end select
    
 end subroutine Flush
  
end module FatesRestartVariableMod
