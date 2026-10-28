module shr_solar_diag_mod
  use shr_kind_mod, only: r8 => SHR_KIND_R8
  use shr_sys_mod, only: shr_sys_abort
  use shr_log_mod, only: shr_log_Unit
  use shr_const_mod, only: shr_const_pi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: shr_solar_diag_configure, shr_solar_diag_schedule
  public :: shr_solar_diag_import, shr_solar_diag_consume
  public :: shr_solar_diag_column, shr_solar_diag_buffer
  public :: shr_solar_diag_radiation
  public :: shr_solar_diag_clock
  logical, public, protected :: shr_solar_diag_enabled = .false.
  logical, public, protected :: shr_solar_diag_requested = .false.
  integer, public, protected :: shr_solar_diag_ymd = 0
  real(r8), public, protected :: shr_solar_diag_tod = -1._r8, shr_solar_diag_dt = 0._r8

  type albedo_record
      ! Rows are visible direct/diffuse, then near-IR direct/diffuse; columns are physics columns.
     real(r8), allocatable :: values(:,:)
     integer :: ymd = 0
     real(r8) :: tod = -1._r8
     logical :: valid = .false.
  end type
  type(albedo_record), allocatable :: imports(:)
  integer :: scheduled_ymd = 0
  real(r8) :: scheduled_tod = -1._r8
  logical :: scheduled_valid = .false.

contains

  subroutine shr_solar_diag_configure(enabled,initialize_only)
    ! Register native restart fields before component initialization, but trace only in cime_run.
    logical, intent(in) :: enabled
    logical, optional, intent(in) :: initialize_only
    shr_solar_diag_requested = enabled
    shr_solar_diag_enabled = enabled
    if (present(initialize_only)) then
       if (initialize_only) shr_solar_diag_enabled = .false.
    endif
    if (.not. enabled .and. allocated(imports)) deallocate(imports)
    scheduled_valid = .false.
  end subroutine

  subroutine shr_solar_diag_clock(ymd,tod,dt)
    integer, intent(in) :: ymd
    real(r8), intent(in) :: tod,dt
    if (.not. shr_solar_diag_enabled) return
    if (.not. ieee_is_finite(tod) .or. .not. ieee_is_finite(dt)) call shr_sys_abort('SOLAR_DIAG: nonfinite clock')
    if (tod < 0._r8 .or. tod >= 86400._r8 .or. dt <= 0._r8) call shr_sys_abort('SOLAR_DIAG: invalid clock')
    shr_solar_diag_ymd = ymd
    shr_solar_diag_tod = tod
    shr_solar_diag_dt = dt
  end subroutine

  subroutine shr_solar_diag_schedule(ymd,tod,valid)
    integer, intent(in) :: ymd
    real(r8), intent(in) :: tod
    logical, intent(in) :: valid
    if (.not. shr_solar_diag_enabled) return
    scheduled_ymd = ymd
    scheduled_tod = tod
    scheduled_valid = valid
  end subroutine

  subroutine shr_solar_diag_import(chunk,values,last_chunk)
    integer, intent(in) :: chunk
    integer, optional, intent(in) :: last_chunk
    real(r8), intent(in) :: values(:,:)
    type(albedo_record), allocatable :: grown(:)
    integer :: upper, saved
    if (.not. shr_solar_diag_enabled) return
    if (chunk < 1 .or. size(values,1) /= 4) call shr_sys_abort('SOLAR_DIAG: invalid albedo import')
     upper = chunk
     if (present(last_chunk)) then
       if (last_chunk < chunk) call shr_sys_abort('SOLAR_DIAG: invalid albedo chunk range')
       upper = last_chunk
     endif
    ! The importer supplies its local upper bound; fallback callers grow geometrically without copying values.
    if (.not. allocated(imports)) allocate(imports(chunk:upper))
     if (chunk < lbound(imports,1) .or. upper > ubound(imports,1)) then
       allocate(grown(min(chunk,lbound(imports,1)):max(upper,ubound(imports,1)+size(imports))))
       do saved = lbound(imports,1), ubound(imports,1)
         call move_alloc(imports(saved)%values,grown(saved)%values)
         grown(saved)%ymd = imports(saved)%ymd
         grown(saved)%tod = imports(saved)%tod
         grown(saved)%valid = imports(saved)%valid
       enddo
       call move_alloc(grown,imports)
    endif
    imports(chunk)%values = values
    imports(chunk)%ymd = scheduled_ymd
    imports(chunk)%tod = scheduled_tod
    imports(chunk)%valid = scheduled_valid
  end subroutine

  subroutine shr_solar_diag_consume(chunk,ymd,tod,values)
    integer, intent(in) :: chunk,ymd
    real(r8), intent(in) :: tod,values(:,:)
    character(384) :: detail
    integer :: location(2)
    if (.not. shr_solar_diag_enabled) return
    if (.not. allocated(imports)) call shr_sys_abort('SOLAR_DIAG: consumption without albedo import')
    if (chunk < lbound(imports,1) .or. chunk > ubound(imports,1)) call shr_sys_abort('SOLAR_DIAG: unknown albedo chunk')
    if (.not. allocated(imports(chunk)%values)) call shr_sys_abort('SOLAR_DIAG: missing albedo chunk')
    if (.not. imports(chunk)%valid) call shr_sys_abort('SOLAR_DIAG: albedo has no scheduled validity')
        if (imports(chunk)%ymd /= ymd .or. imports(chunk)%tod /= tod) then
        write(detail,'(a,i0,a,i8.8,a,f12.4,a,i8.8,a,f12.4)') &
          'SOLAR_DIAG: stale albedo consumption time; chunk=',chunk, &
          '; consumed=',ymd,':',tod,'; scheduled=',imports(chunk)%ymd,':',imports(chunk)%tod
        call shr_sys_abort(trim(detail))
        endif
    if (any(shape(values) /= shape(imports(chunk)%values))) call shr_sys_abort('SOLAR_DIAG: albedo shape differs')
    if (any(.not. ieee_is_finite(values))) call shr_sys_abort('SOLAR_DIAG: nonfinite consumed albedo')
    if (any(values < 0._r8) .or. any(values > 1._r8)) call shr_sys_abort('SOLAR_DIAG: consumed albedo outside [0,1]')
     if (any(values /= imports(chunk)%values)) then
       location = maxloc(abs(values-imports(chunk)%values))
       write(detail,'(a,i0,a,i8.8,a,f12.4,a,i0,a,i0,a,es24.16,a,es24.16)') &
          'SOLAR_DIAG: imported and consumed albedos differ; chunk=',chunk,'; consumed=',ymd,':',tod, &
          '; column=',location(2),'; channel=',location(1), &
          '; imported=',imports(chunk)%values(location(1),location(2)), &
          '; consumed=',values(location(1),location(2))
       call shr_sys_abort(trim(detail))
     endif
    write(shr_log_Unit,'(a,1x,i8.8,1x,f12.4,2(1x,i8))') 'SOLAR_CONSUME_V1',ymd,tod,chunk,size(values,2)
  end subroutine

  subroutine shr_solar_diag_column(chunk,ymd,tod,toa,surface,heating,area)
    integer, intent(in) :: chunk,ymd
    real(r8), intent(in) :: tod,toa(:),surface(:),heating(:)
    real(r8), optional, intent(in) :: area(:)
    real(r8) :: residual(size(toa)),scale(size(toa)),weights(size(toa))
    integer :: column
    if (.not. shr_solar_diag_enabled) return
    if (size(toa) /= size(surface) .or. size(toa) /= size(heating)) &
         call shr_sys_abort('SOLAR_DIAG: column shapes differ')
    if (any(.not. ieee_is_finite(toa)) .or. any(.not. ieee_is_finite(surface)) .or. &
         any(.not. ieee_is_finite(heating))) call shr_sys_abort('SOLAR_DIAG: nonfinite column flux')
    residual = toa - surface - heating
    scale = abs(toa) + abs(surface) + abs(heating)
     if (any(abs(residual) > 1.e-7_r8 + 1.e-10_r8*scale)) then
       column = maxloc(abs(residual)-(1.e-7_r8+1.e-10_r8*scale),dim=1)
       call solar_diag_closure_abort('TOA/heating closure failed',chunk,ymd,tod,column,residual(column),scale(column))
     endif
     weights = 1._r8
     if (present(area)) then
       if (size(area) /= size(toa)) call shr_sys_abort('SOLAR_DIAG: column area shape differs')
       if (any(.not. ieee_is_finite(area))) call shr_sys_abort('SOLAR_DIAG: nonfinite column area')
       if (any(area <= 0._r8)) call shr_sys_abort('SOLAR_DIAG: nonpositive column area')
       weights = area/(4._r8*shr_const_pi)
     endif
     write(shr_log_Unit,'(a,1x,i8.8,1x,f12.4,2(1x,i8),6(1x,es25.16e3))') &
        'SOLAR_COLUMN_V1',ymd,tod,chunk,size(toa),maxval(abs(residual)), &
        sum(toa*weights),sum(surface*weights),sum(heating*weights),maxval(scale),sum(weights)
  end subroutine

    ! toa is MODEL-TOP net SW; top(:,1) is true net TOA SW and top(:,2) is extra-layer absorption [W/m2].
    ! Integrate qrs [W/kg] against pdel [Pa]/gravity [m/s2]. Area is steradians; omitted area uses unit weights.
    ! Held steps reuse saved SOLAR_TOP and heating; only dosw consumes albedo at its actual valid date/TOD.
    subroutine shr_solar_diag_radiation(chunk,ymd,tod,albedo,toa,surface,qrs,pdel,gravity,dosw,area,top)
    integer, intent(in) :: chunk,ymd
    real(r8), intent(in) :: tod,albedo(:,:),toa(:),surface(:),qrs(:,:),pdel(:,:),gravity
    logical, intent(in) :: dosw
    real(r8), optional, intent(in) :: area(:),top(:,:)
    real(r8) :: residual(size(toa)),scale(size(toa))
    integer :: column
    if (.not. shr_solar_diag_enabled) return
    if (.not. ieee_is_finite(gravity)) call shr_sys_abort('SOLAR_DIAG: invalid gravity')
    if (gravity <= 0._r8) call shr_sys_abort('SOLAR_DIAG: nonpositive gravity')
    if (any(shape(qrs) /= shape(pdel))) call shr_sys_abort('SOLAR_DIAG: heating pressure shapes differ')
    if (dosw) call shr_solar_diag_consume(chunk,ymd,tod,albedo)
    if (present(top)) then
       if (size(top,1) /= size(toa) .or. size(top,2) /= 2) call shr_sys_abort('SOLAR_DIAG: upper layer shape differs')
       if (any(.not. ieee_is_finite(top))) call shr_sys_abort('SOLAR_DIAG: nonfinite upper layer')
       residual = top(:,1) - toa - top(:,2)
       scale = abs(top(:,1)) + abs(toa) + abs(top(:,2))
       if (any(abs(residual) > 1.e-7_r8 + 1.e-10_r8*scale)) then
         column = maxloc(abs(residual)-(1.e-7_r8+1.e-10_r8*scale),dim=1)
         call solar_diag_closure_abort('upper layer closure failed',chunk,ymd,tod,column,residual(column),scale(column))
       endif
       call shr_solar_diag_column(chunk,ymd,tod,top(:,1),surface,sum(qrs*pdel,dim=2)/gravity+top(:,2),area)
    else
       call shr_solar_diag_column(chunk,ymd,tod,toa,surface,sum(qrs*pdel,dim=2)/gravity,area)
    endif
  end subroutine

  subroutine solar_diag_closure_abort(message,chunk,ymd,tod,column,residual,scale)
    character(*), intent(in) :: message
    integer, intent(in) :: chunk,ymd,column
    real(r8), intent(in) :: tod,residual,scale
    character(384) :: detail
    write(detail,'(a,a,a,i0,a,i8.8,a,f12.4,a,i0,3(a,es24.16))') &
         'SOLAR_DIAG: ',message,'; chunk=',chunk,'; date=',ymd,':',tod,'; column=',column, &
         '; residual=',residual,'; scale=',scale,'; tolerance=',1.e-7_r8+1.e-10_r8*scale
    call shr_sys_abort(trim(detail))
  end subroutine

  subroutine shr_solar_diag_buffer(before,incoming,delivered,after)
    real(r8), intent(in) :: before(:),incoming(:),delivered(:),after(:)
    real(r8) :: residual(size(before)),scale(size(before))
    if (.not. shr_solar_diag_enabled) return
    if (size(incoming) /= size(before) .or. size(delivered) /= size(before) .or. size(after) /= size(before)) &
         call shr_sys_abort('SOLAR_DIAG: buffer shapes differ')
    if (any(.not. ieee_is_finite(before)) .or. any(.not. ieee_is_finite(incoming)) .or. &
         any(.not. ieee_is_finite(delivered)) .or. any(.not. ieee_is_finite(after))) &
         call shr_sys_abort('SOLAR_DIAG: nonfinite buffer inventory')
    residual = before + incoming - delivered - after
    scale = abs(before) + abs(incoming) + abs(delivered) + abs(after)
    if (any(abs(residual) > 1.e-7_r8 + 1.e-10_r8*scale)) call shr_sys_abort('SOLAR_DIAG: buffer energy closure failed')
  end subroutine

end module shr_solar_diag_mod