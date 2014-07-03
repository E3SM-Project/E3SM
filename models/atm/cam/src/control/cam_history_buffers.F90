module cam_history_buffers
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_history_support, only: max_chars, field_info, hentry
  use abortutils, only : endrun
  use pio, only : var_desc_t
  
  implicit none
!
! dim_index_2d, dim_index_3d: 2-D & 3-D dimension index lower & upper bounds
!
  type dim_index_2d                   ! 2-D dimension index
     sequence
     integer :: beg1, end1            ! lower & upper bounds of 1st dimension
     integer :: beg2, end2            ! lower & upper bounds of 2nd dimension
  end type dim_index_2d
  
  type dim_index_3d                   ! 3-D dimension index
     integer :: beg1, end1            ! lower & upper bounds of 1st dimension
     integer :: beg2, end2            ! lower & upper bounds of 2nd dimension
     integer :: beg3, end3            ! lower & upper bounds of 3rd dimension
  end type dim_index_3d


contains
  subroutine hbuf_accum_inst (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Accumulate instantaneous values of field in 2-D hbuf.
    !          Set accumulation counter to 1.
    !
    !-----------------------------------------------------------------------
    !
    real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
    integer, pointer                 :: nacs(:) ! accumulation counter
    integer, intent(in)              :: idim    ! Longitude dimension of field array
    logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
    real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
    type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
    real(r8), intent(in) :: fillvalue
    !
    ! Local indices
    !
    integer :: ib, ie    ! beginning and ending indices of first dimension
    integer :: jb, je    ! beginning and ending indices of second dimension
    integer :: ieu, jeu  ! number of elements in each dimension
    integer :: i, k      ! loop indices

    logical :: bad       ! flag indicates input field fillvalues not applied consistently
    ! with vertical level

    ib = dimind%beg1
    ie = dimind%end1
    jb = dimind%beg2
    je = dimind%end2

    ieu = ie-ib+1
    jeu = je-jb+1

    do k=1,jeu
       do i=1,ieu
          buf8(i,k) = field(i,k)
       end do
    end do

    if (flag_xyfill) then
       do i=1,ieu
          if (field(i,1) == fillvalue) then
             nacs(i) = 0
          else
             nacs(i) = 1
          end if
       end do
    else
       nacs(1) = 1
    end if

    return
  end subroutine hbuf_accum_inst
  !#######################################################################

  subroutine hbuf_accum_add (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Add the values of field to 2-D hbuf.
    !          Increment accumulation counter by 1.
    !
    !-----------------------------------------------------------------------
    !
    real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
    integer, pointer                 :: nacs(:) ! accumulation counter
    integer, intent(in) :: idim           ! Longitude dimension of field array
    logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
    real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
    type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
    real(r8), intent(in) :: fillvalue
    !
    ! Local indices
    !
    integer :: ib, ie    ! beginning and ending indices of first dimension
    integer :: jb, je    ! beginning and ending indices of second dimension
    integer :: ieu, jeu  ! number of elements in each dimension
    integer :: i,k       ! indices

    ib = dimind%beg1
    ie = dimind%end1
    jb = dimind%beg2
    je = dimind%end2

    ieu = ie-ib+1
    jeu = je-jb+1

    if (flag_xyfill) then
       do k=1,jeu
          do i=1,ieu
             if (field(i,k) /= fillvalue) then
                buf8(i,k) = buf8(i,k) + field(i,k)
             end if
          end do
       end do
       !
       ! Ensure input field has fillvalue defined invariant in the z-direction, then increment nacs
       !
       call check_accum (field, idim, ieu, jeu, fillvalue)
       do i=1,ieu
          if (field(i,1) /= fillvalue) then
             nacs(i) = nacs(i) + 1
          end if
       end do
    else
       do k=1,jeu
          do i=1,ieu
             buf8(i,k) = buf8(i,k) + field(i,k)
          end do
       end do
       nacs(1) = nacs(1) + 1
    end if

    return
  end subroutine hbuf_accum_add

  !#######################################################################

  subroutine hbuf_accum_add00z (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Add the values of field to 2-D hbuf, only of time is 00z.
    !          Increment accumulation counter by 1.
    !
    !-----------------------------------------------------------------------
    !
    use time_manager, only:  get_curr_date

    real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
    integer, pointer                 :: nacs(:) ! accumulation counter
    integer, intent(in) :: idim           ! Longitude dimension of field array
    logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
    real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
    type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
    real(r8), intent(in) :: fillvalue
    !
    ! Local indices
    !
    integer :: ib, ie    ! beginning and ending indices of first dimension
    integer :: jb, je    ! beginning and ending indices of second dimension
    integer :: ieu, jeu  ! number of elements in each dimension
    integer :: i,k       ! indices
    integer :: yr, mon, day, tod

    ! get the time of day, return if not 00z
    call get_curr_date (yr,mon,day,tod)
    if (tod /= 0) return

    ib = dimind%beg1
    ie = dimind%end1
    jb = dimind%beg2
    je = dimind%end2

    ieu = ie-ib+1
    jeu = je-jb+1

    if (flag_xyfill) then
       do k=1,jeu
          do i=1,ieu
             if (field(i,k) /= fillvalue) then
                buf8(i,k) = buf8(i,k) + field(i,k)
             end if
          end do
       end do
       !
       ! Ensure input field has fillvalue defined invariant in the z-direction, then increment nacs
       !
       call check_accum (field, idim, ieu, jeu, fillvalue)
       do i=1,ieu
          if (field(i,1) /= fillvalue) then
             nacs(i) = nacs(i) + 1
          end if
       end do
    else
       do k=1,jeu
          do i=1,ieu
             buf8(i,k) = buf8(i,k) + field(i,k)
          end do
       end do
       nacs(1) = nacs(1) + 1
    end if

    return
  end subroutine hbuf_accum_add00z

  !#######################################################################

  subroutine hbuf_accum_max (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Accumulate the maximum values of field in 2-D hbuf
    !          Set accumulation counter to 1.
    !
    !-----------------------------------------------------------------------
    !
    real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
    integer, pointer                 :: nacs(:) ! accumulation counter
    integer, intent(in) :: idim           ! Longitude dimension of field array
    logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
    real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
    type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
    real(r8), intent(in) :: fillvalue
    !
    ! Local indices
    !
    integer :: ib, ie    ! beginning and ending indices of first dimension
    integer :: jb, je    ! beginning and ending indices of second dimension
    integer :: ieu, jeu  ! number of elements in each dimension
    integer :: i, k

    ib = dimind%beg1
    ie = dimind%end1
    jb = dimind%beg2
    je = dimind%end2

    ieu = ie-ib+1
    jeu = je-jb+1


    if (flag_xyfill) then
       do k=1,jeu
          do i=1,ieu
             if (nacs(i) == 0) then
                buf8(i,k) = -huge (buf8)
             end if
             if (field(i,k) > buf8(i,k) .and. field(i,k) /= fillvalue) then
                buf8(i,k) = field(i,k)
             end if
          end do
       end do
    else
       do k=1,jeu
          do i=1,ieu
             if (nacs(1) == 0) then
                buf8(i,k) = -huge (buf8)
             end if
             if (field(i,k) > buf8(i,k)) then
                buf8(i,k) = field(i,k)
             end if
          end do
       end do
    end if

    if (flag_xyfill) then
       call check_accum (field, idim, ieu, jeu,fillvalue)
       do i=1,ieu
          if (field(i,1) /= fillvalue) then
             nacs(i) = 1
          end if
       end do
    else
       nacs(1) = 1
    end if

    return
  end subroutine hbuf_accum_max

  !#######################################################################

  subroutine hbuf_accum_min (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Accumulate the minimum values of field in 2-D hbuf
    !          Set accumulation counter to 1.
    !
    !-----------------------------------------------------------------------
    !
    real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
    integer, pointer                 :: nacs(:) ! accumulation counter
    integer, intent(in) :: idim           ! Longitude dimension of field array
    logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
    real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
    type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
    real(r8), intent(in) :: fillvalue
    !
    ! Local indices
    !
    integer :: ib, ie    ! beginning and ending indices of first dimension
    integer :: jb, je    ! beginning and ending indices of second dimension
    integer :: ieu, jeu  ! number of elements in each dimension
    integer :: i, k

    ib = dimind%beg1
    ie = dimind%end1
    jb = dimind%beg2
    je = dimind%end2

    ieu = ie-ib+1
    jeu = je-jb+1


    if (flag_xyfill) then
       do k=1,jeu
          do i=1,ieu
             if (nacs(i) == 0) then
                buf8(i,k) = +huge (buf8)
             end if
             if (field(i,k) < buf8(i,k) .and. field(i,k) /= fillvalue) then
                buf8(i,k) = field(i,k)
             end if
          end do
       end do
    else
       do k=1,jeu
          do i=1,ieu
             if (nacs(1) == 0) then
                buf8(i,k) = +huge (buf8)
             end if
             if (field(i,k) < buf8(i,k)) then
                buf8(i,k) = field(i,k)
             end if
          end do
       end do
    end if

    if (flag_xyfill) then
       call check_accum (field, idim, ieu, jeu, fillvalue)
       do i=1,ieu
          if (field(i,1) /= fillvalue) then
             nacs(i) = 1
          end if
       end do
    else
       nacs(1) = 1
    end if

    return
  end subroutine hbuf_accum_min

  subroutine hbuf_accum_addlcltime (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue, c , decomp_type,&
       lcltod_start, lcltod_stop)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Add the values of field to 2-D hbuf, only if the local time
    !          is in the range specified.
    !          Increment accumulation counter by 1.
    !
    !-----------------------------------------------------------------------
    !
    use time_manager, only:  get_curr_date
    use phys_grid,    only:  get_rlon_all_p
    use physconst,    only:  pi
    use phys_grid,     only: get_ncols_p, pcols
    use cam_pio_utils, only: phys_decomp
    use dyn_grid,      only: get_horiz_grid_dim_d, get_horiz_grid_d, dyn_grid_get_elem_coords

    type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
    real(r8), pointer    :: buf8(:,:)    ! 2-D history buffer
    integer,  pointer    :: nacs(:) ! accumulation counter
    integer,  intent(in) :: idim           ! Longitude dimension of field array
    logical,  intent(in) :: flag_xyfill ! non-applicable xy points flagged with fillvalue
    real(r8), intent(in) :: field(idim,*)   ! real*8 array
    integer,  intent(in) :: c              ! chunk (physics) or latitude (dynamics) index

    integer, intent(in)  :: decomp_type, lcltod_start, lcltod_stop
    real(r8), intent(in) :: fillvalue
    
    !
    ! Local indices
    !
    integer  :: ib, ie    ! beginning and ending indices of first dimension
    integer  :: jb, je    ! beginning and ending indices of second dimension
    integer  :: ieu, jeu  ! number of elements in each dimension
    integer  :: i,k       ! indices
    integer  :: yr, mon, day, tod
    integer  :: ncols, hdim1_d, hdim2_d

    integer,  allocatable :: lcltod(:)     ! local time of day (secs)
    logical,  allocatable :: inavg(:)      ! is the column in the desired local time range?
    real(r8), allocatable :: rlon(:)       ! column longitude (radians)
    integer,  allocatable :: cdex(:)       ! global column index 

    ib = dimind%beg1
    ie = dimind%end1
    jb = dimind%beg2
    je = dimind%end2

    ieu = ie-ib+1
    jeu = je-jb+1

    allocate( inavg(1:ieu) , lcltod(1:ieu) )
    lcltod(:) = 0

    !
    ! Get the time of day and longitude and compute the local time.
    !
    call get_curr_date (yr,mon,day,tod)      

    if ( decomp_type == phys_decomp ) then 

       ncols = get_ncols_p(c)
       ie = ncols
       ieu = ncols
       allocate( rlon(ncols) )
       call get_rlon_all_p(c, ncols, rlon)
       lcltod(1:ncols) = mod((tod) + (nint(86400_r8 * rlon(1:ncols) / 2._r8 / pi)), 86400)

    else 

       ncols = ieu
       allocate(rlon(ncols),cdex(ncols))
       call dyn_grid_get_elem_coords( c, rlon=rlon, cdex=cdex )
       lcltod(:) = -999999
       where( cdex(:)>0 ) lcltod(1:ieu) = mod((tod) + (nint(86400_r8 * rlon(1:ncols) / 2._r8 / pi)), 86400)

    endif



    !
    ! Set a flag to indicate that the column is in the requested local time range.      
    ! If lcltod_stop is less than lcltod_stop, then the time is wrapping around 24 hours.
    !
    inavg(:)  = .false.

    if (lcltod_stop < lcltod_start) then
       ! the ".and.(lcltod(:)>0" condition was added to exclude the undifined (-999999) columns
       where((lcltod(:) >= lcltod_start) .or. ((lcltod(:) < lcltod_stop).and.(lcltod(:)>0))) inavg(:) = .true.
    else if (lcltod_stop == lcltod_start) then
       where(lcltod(ib:ie) == lcltod_start) inavg(1:ieu) = .true.
    else
       where((lcltod(:) >= lcltod_start) .and. (lcltod(:) < lcltod_stop)) inavg(:) = .true.
    end if

    if (flag_xyfill) then
       do k=1,jeu
          do i=1,ieu
             if (inavg(i) .and. (field(i,k) /= fillvalue)) then
                buf8(i,k) = buf8(i,k) + field(i,k)
             end if
          end do
       end do
       !
       ! Ensure input field has fillvalue defined invariant in the z-direction, then increment nacs
       !
       call check_accum (field, idim, ieu, jeu, fillvalue)
       do i=1,ieu
          if (inavg(i) .and. (field(i,1) /= fillvalue)) then
             nacs(i) = nacs(i) + 1
          end if
       end do
    else
       do k=1,jeu
          do i=1,ieu
             if (inavg(i)) then
                buf8(i,k) = buf8(i,k) + field(i,k)
             end if
          end do
       end do
       !
       ! NOTE: Because of the local time, some columns in the chunk may not be used in the
       ! average, so nacs need to be dimension with the number of columns unlike the other 
       ! averaging techniques.
       !
       do i=1,ieu
          if (inavg(i)) nacs(i) = nacs(i) + 1
       end do
    end if

    deallocate( inavg , rlon, lcltod )
    if (allocated(cdex)) deallocate(cdex)

    return
  end subroutine hbuf_accum_addlcltime

!#######################################################################


  !#######################################################################

  subroutine check_accum (field, idim, ieu, jeu, fillvalue)
    !
    integer, intent(in)  :: idim
    real(r8), intent(in) :: field(idim,*)   ! real*8 array
    integer, intent(in)  :: ieu,jeu         ! loop ranges

    logical :: bad
    integer :: i,k
    real(r8), intent(in) :: fillvalue
    !
    ! For multilevel fields ensure that all levels have fillvalue applied consistently
    !
    bad = .false.
    do k=2,jeu
       do i=1,ieu
          if (field(i,1) == fillvalue .and. field(i,k) /= fillvalue .or. &
               field(i,1) /= fillvalue .and. field(i,k) == fillvalue) then
             bad = .true.
          end if
       end do
    end do

    if (bad) then
       call endrun ('CHECK_ACCUM: inconsistent level values')
    end if

    return
  end subroutine check_accum

end module cam_history_buffers
