module accumulMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains generic subroutines that can be used to
  ! define, accumulate and extract  user-specified fields over
  ! user-defined intervals. Each interval  and accumulation type is
  ! unique to each field processed.
  ! Subroutine [init_accumulator] defines the values of the accumulated
  ! field data structure. Subroutine [update_accum_field] does
  ! the actual accumulation for a given field.
  ! Four types of accumulations are possible:
  ! - Average over time interval. Time average fields are only
  !   valid at the end of the averaging interval.
  ! - Running mean over time interval. Running means are valid once the
  !   length of the simulation exceeds the
  ! - Running accumulation over time interval. Accumulated fields are
  !   continuously accumulated. The trigger value "-99999." resets
  !   the accumulation to zero.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_sys_mod , only: shr_sys_abort
  use clm_varctl  , only: iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: accumulRest          ! Write/read restart of accumulation fields
  public :: init_accum_field     ! Initialize an accumulator field
  public :: print_accum_fields   ! Print info about accumulator fields
  public :: extract_accum_field  ! Extracts the current value of an accumulator field
  public :: update_accum_field   ! Update the current value of an accumulator field

  interface extract_accum_field
     module procedure extract_accum_field_sl ! Extract current val of single-level accumulator field
     module procedure extract_accum_field_ml ! Extract current val of multi-level accumulator field
  end interface
  interface update_accum_field               ! Updates the current value of an accumulator field
     module procedure update_accum_field_sl  ! Update single-level accumulator field
     module procedure update_accum_field_ml  ! Update multi-level accumulator field
  end interface
  private
  !
  ! !PRIVATE TYPES:
  type accum_field
     character(len=  8) :: name     !field name
     character(len=128) :: desc     !field description
     character(len=  8) :: units    !field units
     character(len=  8) :: acctype  !accumulation type: ["timeavg","runmean","runaccum"]
     character(len=  8) :: type1d   !subgrid type: ["gridcell","landunit","column" or "pft"]
     character(len=  8) :: type2d   !type2d ('','levsoi','numrad',..etc. )
     integer            :: beg1d    !subgrid type beginning index
     integer            :: end1d    !subgrid type ending index
     integer            :: num1d    !total subgrid points
     integer            :: numlev   !number of vertical levels in field
     real(r8)           :: initval  !initial value of accumulated field
     real(r8), pointer  :: val(:,:) !accumulated field
     integer            :: period   !field accumulation period (in model time steps)
  end type accum_field

  real(r8), parameter, public :: accumResetVal = -99999._r8 ! used to do an annual reset ( put in for bug 1858)

  integer, parameter :: max_accum = 100    !maximum number of accumulated fields
  type (accum_field) :: accum(max_accum)   !array accumulated fields
  integer :: naccflds = 0                  !accumulator field counter

  integer :: iflag_interp = 1
  integer :: iflag_copy   = 2
  integer :: iflag_skip   = 3
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine init_accum_field (name, units, desc, &
       accum_type, accum_period, numlev, subgrid_type, init_value, type2d)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation fields. This subroutine sets:
    ! o name  of accumulated field
    ! o units of accumulated field
    ! o accumulation type of accumulated field
    ! o description of accumulated fields: accdes
    ! o accumulation period for accumulated field (in iterations)
    ! o initial value of accumulated field
    !
    ! !USES:
    use shr_const_mod, only: SHR_CONST_CDAY
    use clm_time_manager, only : get_step_size
    use decompMod, only : get_proc_bounds, get_proc_global
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)           :: name         !field name
    character(len=*), intent(in)           :: units        !field units
    character(len=*), intent(in)           :: desc         !field description
    character(len=*), intent(in)           :: accum_type   !field type: tavg, runm, runa, ins
    integer , intent(in)                   :: accum_period !field accumulation period
    character(len=*), intent(in)           :: subgrid_type !["gridcell","landunit","column" or "pft"]
    integer , intent(in)                   :: numlev       !number of vertical levels
    real(r8), intent(in)                   :: init_value   !field initial or reset value
    character(len=*), intent(in), optional :: type2d       !level type (optional) - needed if numlev > 1
    !
    ! !LOCAL VARIABLES:
    integer :: nf           ! field index
    integer :: beg1d,end1d  ! beggining and end subgrid indices
    integer :: num1d        ! total number subgrid indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: begCohort, endCohort   ! per-proc beg end cohort indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    integer :: numCohort    ! total number of cohorts across all processors
    !------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp, &
         begCohort, endCohort )
    call get_proc_global(numg, numl, numc, nump, numCohort)

    ! update field index
    ! Consistency check that number of accumulated does not exceed maximum.

    naccflds = naccflds + 1
    if (naccflds > max_accum) then
       write(iulog,*) 'ACCUMULINIT error: user-defined accumulation fields ', &
            'equal to ',naccflds,' exceeds max_accum'
       call shr_sys_abort()
    end if
    nf = naccflds

    ! Note accumulation period must be converted from days
    ! to number of iterations

    accum(nf)%name    = trim(name)
    accum(nf)%units   = trim(units)
    accum(nf)%desc    = trim(desc)
    accum(nf)%acctype = trim(accum_type)
    accum(nf)%initval = init_value
    accum(nf)%period  = accum_period
    if (accum(nf)%period < 0) then
       accum(nf)%period = -accum(nf)%period * nint(SHR_CONST_CDAY) / get_step_size()
    end if

    select case (trim(subgrid_type))
    case ('gridcell')
       beg1d = begg
       end1d = endg
       num1d = numg
    case ('landunit')
       beg1d = begl
       end1d = endl
       num1d = numl
    case ('column')
       beg1d = begc
       end1d = endc
       num1d = numc
    case ('pft')
       beg1d = begp
       end1d = endp
       num1d = nump
    case default
       write(iulog,*)'ACCUMULINIT: unknown subgrid type ',subgrid_type
       call shr_sys_abort ()
    end select

    accum(nf)%type1d = trim(subgrid_type)
    accum(nf)%beg1d = beg1d
    accum(nf)%end1d = end1d
    accum(nf)%num1d = num1d
    accum(nf)%numlev = numlev

    if (present(type2d)) then
       accum(nf)%type2d = type2d
    else
       accum(nf)%type2d = ' '
    end if
    
    ! Allocate and initialize accumulation field

    allocate(accum(nf)%val(beg1d:end1d,numlev))
    accum(nf)%val(beg1d:end1d,numlev) = init_value

  end subroutine init_accum_field

  !------------------------------------------------------------------------
  subroutine print_accum_fields()
    !
    ! !DESCRIPTION:
    ! Diagnostic printout of accumulated fields
    !
    ! !USES:
    use spmdMod, only : masterproc
    !
    ! !ARGUMENTS:
    implicit none
    !
    integer :: i,nf   !indices
    !------------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*)
       write(iulog,*) 'Initializing variables for time accumulation .....'
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Accumulated fields'
       write(iulog,1002)
       write(iulog,'(72a1)') ("_",i=1,71)
       do nf = 1, naccflds
          if (accum(nf)%period /= huge(1)) then
             write(iulog,1003) nf,accum(nf)%name,accum(nf)%units,&
                  accum(nf)%acctype, accum(nf)%period, accum(nf)%initval, &
                  accum(nf)%desc
          else
             write(iulog,1004) nf,accum(nf)%name,accum(nf)%units,&
                  accum(nf)%acctype, accum(nf)%initval, accum(nf)%desc
          endif
       end do
       write(iulog,'(72a1)') ("_",i=1,71)
       write(iulog,*)
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Successfully initialized variables for accumulation'
       write(iulog,*)
    endif

1002 format(' No',' Name    ',' Units   ',' Type    ','Period',' Inival',' Description')
1003 format((1x,i2),(1x,a8),(1x,a8),(1x,a8), (1x,i5),(1x,f4.0),(1x,a40))
1004 format((1x,i2),(1x,a8),(1x,a8),(1x,a8),'  N.A.',(1x,f4.0),(1x,a40))

  end subroutine print_accum_fields

  !------------------------------------------------------------------------
  subroutine extract_accum_field_sl (name, field, nstep)
    !
    ! !DESCRIPTION:
    ! Extract single-level accumulated field.
    ! This routine extracts the field values from the multi-level
    ! accumulation field. It extracts the current value except if
    ! the field type is a time average. In this case, an absurd value
    ! is assigned to  indicate the time average is not yet valid.
    !
    ! !USES:
    use clm_varcon, only : spval, ispval
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name     !field name
    real(r8), pointer, dimension(:) :: field !field values for current time step
    integer , intent(in) :: nstep            !timestep index
    !
    ! !LOCAL VARIABLES:
    integer :: i,k,nf        !indices
    integer :: beg,end         !subgrid beginning,ending indices
    !------------------------------------------------------------------------

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1, naccflds
       if (name == accum(i)%name) nf = i
    end do
    if (nf == 0) then
       write(iulog,*) 'EXTRACT_ACCUM_FIELD_SL error: field name ',name,' not found'
       call shr_sys_abort
    endif

    ! error check

    beg = accum(nf)%beg1d
    end = accum(nf)%end1d
    if (size(field,dim=1) /= end-beg+1) then
       write(iulog,*)'ERROR in extract_accum_field for field ',accum(nf)%name
       write(iulog,*)'size of first dimension of field is ',&
            size(field,dim=1),' and should be ',end-beg+1
       call shr_sys_abort
    endif

    ! extract field

    if (accum(nf)%acctype == 'timeavg' .and. &
         mod(nstep,accum(nf)%period) /= 0) then
       do k = beg,end
          field(k) = spval  !assign absurd value when avg not ready
       end do
    else
       do k = beg,end
          field(k) = accum(nf)%val(k,1)
       end do
    end if

  end subroutine extract_accum_field_sl

  !------------------------------------------------------------------------
  subroutine extract_accum_field_ml (name, field, nstep)
    !
    ! !DESCRIPTION:
    ! Extract mutli-level accumulated field.
    ! This routine extracts the field values from the multi-level
    ! accumulation field. It extracts the current value except if
    ! the field type is a time average. In this case, an absurd value
    ! is assigned to  indicate the time average is not yet valid.
    !
    ! !USES:
    use clm_varcon, only : spval
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name       !field name
    real(r8), pointer, dimension(:,:) :: field !field values for current time step
    integer, intent(in) :: nstep               !timestep index
    !
    ! !LOCAL VARIABLES:
    integer :: i,j,k,nf        !indices
    integer :: beg,end         !subgrid beginning,ending indices
    integer :: numlev          !number of vertical levels
    !------------------------------------------------------------------------

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1, naccflds
       if (name == accum(i)%name) nf = i
    end do
    if (nf == 0) then
       write(iulog,*) 'EXTRACT_ACCUM_FIELD_ML error: field name ',name,' not found'
       call shr_sys_abort
    endif

    ! error check

    numlev = accum(nf)%numlev
    beg = accum(nf)%beg1d
    end = accum(nf)%end1d
    if (size(field,dim=1) /= end-beg+1) then
       write(iulog,*)'ERROR in extract_accum_field for field ',accum(nf)%name
       write(iulog,*)'size of first dimension of field is ',&
            size(field,dim=1),' and should be ',end-beg+1
       call shr_sys_abort
    else if (size(field,dim=2) /= numlev) then
       write(iulog,*)'ERROR in extract_accum_field for field ',accum(nf)%name
       write(iulog,*)'size of second dimension of field iis ',&
            size(field,dim=2),' and should be ',numlev
       call shr_sys_abort
    endif

    !extract field

    if (accum(nf)%acctype == 'timeavg' .and. &
         mod(nstep,accum(nf)%period) /= 0) then
       do j = 1,numlev
          do k = beg,end
             field(k,j) = spval  !assign absurd value when avg not ready
          end do
       end do
    else
       do j = 1,numlev
          do k = beg,end
             field(k,j) = accum(nf)%val(k,j)
          end do
       end do
    end if

  end subroutine extract_accum_field_ml

  !------------------------------------------------------------------------
  subroutine update_accum_field_sl (name, field, nstep)
    !
    ! !DESCRIPTION:
    ! Accumulate single level field over specified time interval.
    ! The appropriate field is accumulated in the array [accval].
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name     !field name
    real(r8), pointer, dimension(:) :: field !field values for current time step
    integer , intent(in) :: nstep            !time step index
    !
    ! !LOCAL VARIABLES:
    integer :: i,k,nf              !indices
    integer :: accper              !temporary accumulation period
    integer :: beg,end             !subgrid beginning,ending indices
    !------------------------------------------------------------------------

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1, naccflds
       if (name == accum(i)%name) nf = i
    end do
    if (nf == 0) then
       write(iulog,*) 'UPDATE_ACCUM_FIELD_SL error: field name ',name,' not found'
       call shr_sys_abort
    endif

    ! error check

    beg = accum(nf)%beg1d
    end = accum(nf)%end1d
    if (size(field,dim=1) /= end-beg+1) then
       write(iulog,*)'ERROR in UPDATE_ACCUM_FIELD_SL for field ',accum(nf)%name
       write(iulog,*)'size of first dimension of field is ',size(field,dim=1),&
            ' and should be ',end-beg+1
       call shr_sys_abort
    endif

    ! accumulate field

    if (accum(nf)%acctype /= 'timeavg' .AND. &
        accum(nf)%acctype /= 'runmean' .AND. &
        accum(nf)%acctype /= 'runaccum') then
       write(iulog,*) 'UPDATE_ACCUM_FIELD_SL error: incorrect accumulation type'
       write(iulog,*) ' was specified for field ',name
       write(iulog,*)' accumulation type specified is ',accum(nf)%acctype
       write(iulog,*)' only [timeavg, runmean, runaccum] are currently acceptable'
       call shr_sys_abort()
    end if


    ! reset accumulated field value if necessary and  update
    ! accumulation field
    ! running mean never reset

    if (accum(nf)%acctype == 'timeavg') then

       !time average field reset every accumulation period
       !normalize at end of accumulation period

       if ((mod(nstep,accum(nf)%period) == 1 .or. accum(nf)%period == 1) .and. (nstep /= 0))then
          accum(nf)%val(beg:end,1) = 0._r8
       end if
       accum(nf)%val(beg:end,1) =  accum(nf)%val(beg:end,1) + field(beg:end)
       if (mod(nstep,accum(nf)%period) == 0) then
          accum(nf)%val(beg:end,1) = accum(nf)%val(beg:end,1) / accum(nf)%period
       endif

    else if (accum(nf)%acctype == 'runmean') then

       !running mean - reset accumulation period until greater than nstep

       accper = min (nstep,accum(nf)%period)
       accum(nf)%val(beg:end,1) = ((accper-1)*accum(nf)%val(beg:end,1) + field(beg:end)) / accper

    else if (accum(nf)%acctype == 'runaccum') then

       !running accumulation field reset at trigger -99999

       do k = beg,end
          if (nint(field(k)) == -99999) then
             accum(nf)%val(k,1) = 0._r8
          end if
       end do
       accum(nf)%val(beg:end,1) = min(max(accum(nf)%val(beg:end,1) + field(beg:end), 0._r8), 99999._r8)

    end if

  end subroutine update_accum_field_sl

  !------------------------------------------------------------------------
  subroutine update_accum_field_ml (name, field, nstep)
    !
    ! !DESCRIPTION:
    ! Accumulate multi level field over specified time interval.
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name       !field name
    real(r8), pointer, dimension(:,:) :: field !field values for current time step
    integer , intent(in) :: nstep              !time step index
    !
    ! !LOCAL VARIABLES:
    integer :: i,j,k,nf            !indices
    integer :: accper              !temporary accumulation period
    integer :: beg,end             !subgrid beginning,ending indices
    integer :: numlev              !number of vertical levels
    !------------------------------------------------------------------------

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1, naccflds
       if (name == accum(i)%name) nf = i
    end do
    if (nf == 0) then
       write(iulog,*) 'UPDATE_ACCUM_FIELD_ML error: field name ',name,' not found'
       call shr_sys_abort
    endif

    ! error check

    numlev = accum(nf)%numlev
    beg = accum(nf)%beg1d
    end = accum(nf)%end1d
    if (size(field,dim=1) /= end-beg+1) then
       write(iulog,*)'ERROR in UPDATE_ACCUM_FIELD_ML for field ',accum(nf)%name
       write(iulog,*)'size of first dimension of field is ',size(field,dim=1),&
            ' and should be ',end-beg+1
       call shr_sys_abort
    else if (size(field,dim=2) /= numlev) then
       write(iulog,*)'ERROR in UPDATE_ACCUM_FIELD_ML for field ',accum(nf)%name
       write(iulog,*)'size of second dimension of field is ',size(field,dim=2),&
            ' and should be ',numlev
       call shr_sys_abort
    endif

    ! accumulate field

    if (accum(nf)%acctype /= 'timeavg' .AND. &
        accum(nf)%acctype /= 'runmean' .AND. &
        accum(nf)%acctype /= 'runaccum') then
       write(iulog,*) 'UPDATE_ACCUM_FIELD_ML error: incorrect accumulation type'
       write(iulog,*) ' was specified for field ',name
       write(iulog,*)' accumulation type specified is ',accum(nf)%acctype
       write(iulog,*)' only [timeavg, runmean, runaccum] are currently acceptable'
       call shr_sys_abort()
    end if

    ! accumulate field

    ! reset accumulated field value if necessary and  update
    ! accumulation field
    ! running mean never reset

    if (accum(nf)%acctype == 'timeavg') then

       !time average field reset every accumulation period
       !normalize at end of accumulation period

       if ((mod(nstep,accum(nf)%period) == 1 .or. accum(nf)%period == 1) .and. (nstep /= 0))then
          accum(nf)%val(beg:end,1:numlev) = 0._r8
       endif
       accum(nf)%val(beg:end,1:numlev) =  accum(nf)%val(beg:end,1:numlev) + field(beg:end,1:numlev)
       if (mod(nstep,accum(nf)%period) == 0) then
          accum(nf)%val(beg:end,1:numlev) = accum(nf)%val(beg:end,1:numlev) / accum(nf)%period
       endif

    else if (accum(nf)%acctype == 'runmean') then

       !running mean - reset accumulation period until greater than nstep

       accper = min (nstep,accum(nf)%period)
       accum(nf)%val(beg:end,1:numlev) = &
            ((accper-1)*accum(nf)%val(beg:end,1:numlev) + field(beg:end,1:numlev)) / accper

    else if (accum(nf)%acctype == 'runaccum') then

       !running accumulation field reset at trigger -99999

       do j = 1,numlev
          do k = beg,end
             if (nint(field(k,j)) == -99999) then
                accum(nf)%val(k,j) = 0._r8
             end if
          end do
       end do
       accum(nf)%val(beg:end,1:numlev) = &
            min(max(accum(nf)%val(beg:end,1:numlev) + field(beg:end,1:numlev), 0._r8), 99999._r8)

    end if

  end subroutine update_accum_field_ml

  !------------------------------------------------------------------------
  subroutine accumulRest( ncid, flag )
    !
    ! !DESCRIPTION:
    ! Read/write accumulation restart data
    !
    ! !USES:
    use clm_time_manager, only : is_restart
    use clm_varcon      , only : ispval
    use ncdio_pio
    use pio
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid   !netcdf unit
    character(len=*) , intent(in) :: flag   !'define','read', or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: nf,k,j                         ! indices
    integer :: beg1d, end1d                   ! buffer bounds
    integer :: ier                            ! error status
    logical :: readvar                        ! determine if variable is on initial file
    real(r8), pointer  :: rbuf1d(:)           ! temporary 1d buffer
    type(var_desc_t)   :: vardesc             ! local vardesc
    character(len=128) :: varname             ! temporary
    character(len= 32) :: subname='AccumRest' ! subroutine name
    !------------------------------------------------------------------------

    do nf = 1,naccflds

       ! Note = below need to allocate rbuf for single level variables, since
       ! accum(nf)%val is always 2d

       varname = trim(accum(nf)%name) // '_VALUE'
       if (flag == 'define') then
          if (accum(nf)%numlev == 1) then
             call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_double,  &
                  dim1name=accum(nf)%type1d, &
                  long_name=accum(nf)%desc, units=accum(nf)%units) 
             ier = PIO_inq_varid(ncid, trim(varname), vardesc)
             ier = PIO_put_att(ncid, vardesc%varid, 'interpinic_flag', iflag_interp)
          else
             call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_double,  &
                  dim1name=accum(nf)%type1d, dim2name=accum(nf)%type2d, &
                  long_name=accum(nf)%desc, units=accum(nf)%units) 
          end if
       else if (flag == 'read' .or. flag == 'write') then
          if (accum(nf)%numlev == 1) then
             beg1d = accum(nf)%beg1d
             end1d = accum(nf)%end1d
             allocate(rbuf1d(beg1d:end1d))
             if (flag == 'write') then
                rbuf1d(beg1d:end1d) = accum(nf)%val(beg1d:end1d,1) 
             end if
             call ncd_io(varname=varname, data=rbuf1d, &
                  dim1name=accum(nf)%type1d, ncid=ncid, flag=flag, readvar=readvar)
             if (flag == 'read' .and. readvar) then
                accum(nf)%val(beg1d:end1d,1) = rbuf1d(beg1d:end1d)
             end if
             deallocate(rbuf1d)
          else
             call ncd_io(varname=varname, data=accum(nf)%val, &
                  dim1name=accum(nf)%type1d, ncid=ncid, flag=flag, readvar=readvar)
          end if
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call shr_sys_abort()
          end if
       end if

       varname = trim(accum(nf)%name) // '_PERIOD'
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
               long_name='', units='time steps', imissing_value=ispval, &
               ifill_value=huge(1))
          ier = PIO_inq_varid(ncid, trim(varname), vardesc)
          ier = PIO_put_att(ncid, vardesc%varid, 'interpinic_flag', iflag_copy)
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname=varname, data=accum(nf)%period, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call shr_sys_abort()
          end if
       end if

    end do

  end subroutine accumulRest

end module accumulMod
