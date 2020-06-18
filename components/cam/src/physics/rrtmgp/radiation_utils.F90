module radiation_utils

   use shr_kind_mod, only: r8=>shr_kind_r8, cl=>shr_kind_cl
   use assertions, only: assert

   implicit none
   private

   public :: compress_day_columns, expand_day_columns, &
             calculate_heating_rate, clip_values, &
             handle_error

   ! Interface blocks for overloaded procedures
   interface compress_day_columns
      module procedure compress_day_columns_1d, compress_day_columns_2d
   end interface compress_day_columns

   interface expand_day_columns
      module procedure expand_day_columns_1d, expand_day_columns_2d
   end interface expand_day_columns

   interface clip_values
      module procedure clip_values_1d, clip_values_2d, clip_values_3d
   end interface clip_values

   ! Max length for character strings
   integer, parameter :: max_char_len = 512

   ! Name of this module for error messages
   character(len=*), parameter :: module_name = 'radiation_utils'

contains

   !-------------------------------------------------------------------------------
   subroutine compress_day_columns_1d(xcol, xday, day_indices)

      real(r8), intent(in) :: xcol(:)
      real(r8), intent(inout) :: xday(:)
      integer, intent(in) :: day_indices(:)
      integer :: icol, iday
      character(len=32) :: subname = 'compress_day_columns'

      ! Check dimensions
      call assert(size(xday, 1) == size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xcol, 1) >= size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      
      ! Do mapping
      do iday = 1,count(day_indices>0)
         icol = day_indices(iday)
         xday(iday) = xcol(icol)
      end do

   end subroutine
   !-------------------------------------------------------------------------------
   subroutine compress_day_columns_2d(xcol, xday, day_indices)

      real(r8), intent(in) :: xcol(:,:)
      real(r8), intent(inout) :: xday(:,:)
      integer, intent(in) :: day_indices(:)
      integer :: icol, iday
      character(len=32) :: subname = 'compress_day_columns'

      ! Check dimensions
      call assert(size(xday, 1) == size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xcol, 1) >= size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xday, 2) == size(xcol, 2), trim(subname) // ': inconsistent sizes')

      ! Do mapping
      do iday = 1,count(day_indices>0)
         icol = day_indices(iday)
         xday(iday,:) = xcol(icol,:)
      end do

   end subroutine
   !-------------------------------------------------------------------------------
   subroutine expand_day_columns_1d(xday, xcol, day_indices)

      real(r8), intent(in) :: xday(:)
      real(r8), intent(inout) :: xcol(:)
      integer, intent(in) :: day_indices(:)
      integer :: icol, iday
      character(len=32) :: subname = 'expand_day_columns_1d'

      ! Check dimensions
      call assert(size(xday, 1) == size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xcol, 1) >= size(day_indices, 1), trim(subname) // ': inconsistent sizes')

      ! We need to reset to zero because we only populate the daytime columns
      xcol(:) = 0._r8

      ! Do mapping
      do iday = 1,count(day_indices>0)
         icol = day_indices(iday)
         xcol(icol) = xday(iday)
      end do
   end subroutine expand_day_columns_1d
   !-------------------------------------------------------------------------------
   subroutine expand_day_columns_2d(xday, xcol, day_indices)

      real(r8), intent(in) :: xday(:,:)
      real(r8), intent(inout) :: xcol(:,:)
      integer, intent(in) :: day_indices(:)
      integer :: icol, iday
      character(len=32) :: subname = 'expand_day_columns_2d'

      ! Check dimensions
      call assert(size(xday, 1) == size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xcol, 1) >= size(day_indices, 1), trim(subname) // ': inconsistent sizes')
      call assert(size(xday, 2) == size(xcol, 2), trim(subname) // ': inconsistent sizes')

      ! We need to reset to zero because we only populate the daytime columns
      xcol(:,:) = 0._r8

      ! Do mapping
      do iday = 1,count(day_indices>0)
         icol = day_indices(iday)
         xcol(icol,:) = xday(iday,:)
      end do

   end subroutine expand_day_columns_2d

   !-------------------------------------------------------------------------------

   subroutine calculate_heating_rate(flux_up, flux_dn, pint, heating_rate)

      use physconst, only: gravit

      ! Inputs
      real(r8), intent(in), dimension(:,:) :: flux_up, flux_dn, pint

      ! Output heating rate; same size as pint with one fewer vertical level
      real(r8), intent(out) :: heating_rate(:,:)

      ! Loop indices
      integer :: icol, ilev

      ! Everyone needs a name
      character(len=32) :: subname = 'calculate_heating_rate'

      ! Get dimension sizes and make sure arrays conform
      call assert(size(pint,1) == size(flux_up,1), subname // ': sizes do not conform.')
      call assert(size(pint,2) == size(flux_up,2), subname // ': sizes do not conform.')
      call assert(size(heating_rate,1) == size(flux_up,1), subname // ': sizes do not conform.')
      call assert(size(heating_rate,2) == size(flux_up,2)-1, subname // ': sizes do not conform.')

      ! Loop over levels and calculate heating rates; note that the fluxes *should*
      ! be defined at interfaces, so the loop ktop,kbot and grabbing the current
      ! and next value of k should be safe. ktop should be the top interface, and
      ! kbot + 1 should be the bottom interface.
      !
      ! NOTE: to get heating rate in K/day, normally we would use:
      !
      !     H = dF / dp * g * (sec/day) * (1e-5) / (cpair)
      !
      ! Here we just use
      !
      !     H = dF / dp * g
      !
      ! Why? Something to do with convenience with applying the fluxes to the
      ! heating tendency?
      do ilev = 1,size(pint,2)-1
         do icol = 1,size(pint,1)
            heating_rate(icol,ilev) = ( &
               flux_up(icol,ilev+1) - flux_up(icol,ilev) - &
               flux_dn(icol,ilev+1) + flux_dn(icol,ilev) &
            ) * gravit / (pint(icol,ilev+1) - pint(icol,ilev))
         end do
      end do

   end subroutine calculate_heating_rate

   !----------------------------------------------------------------------------

   ! Routines to clip array values if they are outside of an expected range,
   ! defined by min_x and max_x. Allow passing a varname argument, which will be
   ! appended to warning message and tolerance argument that sets a maximum
   ! tolerance above which a warning or error will be thrown. That is, these 
   ! routines will *always* clip values outside expected range and allow the 
   ! simulation to continue, but a warning will be issued if the values fall 
   ! outside the expected range plus the tolerance. In other words, no warning 
   ! is issued when clipping values outside the valid range plus/minus the 
   ! tolerance.
   function clip_values_1d(x, min_x, max_x, varname, tolerance) result(error_message)
      real(r8), intent(inout) :: x(:)
      real(r8), intent(in) :: min_x
      real(r8), intent(in) :: max_x
      character(len=*), intent(in) :: varname
      real(r8), intent(in), optional :: tolerance
      real(r8) :: tolerance_local
      character(max_char_len) :: error_message

      error_message = ''
      tolerance_local = 0._r8
      if (present(tolerance)) then
         tolerance_local = tolerance
      end if

      ! look for values less than threshold
      if (any(x < min_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x < (min_x - tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x < (min_x - tolerance_local)), ' values below threshold ', &
                     '; min = ', minval(x)
         end if
         ! Clip values
         where (x < min_x)
            x = min_x
         endwhere
      end if

      ! Look for values greater than threshold
      if (any(x > max_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x > (max_x + tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x > (max_x + tolerance_local)), ' values above threshold ', &
                     '; max = ', maxval(x)
         end if
         ! Clip values
         where (x > max_x)
            x = max_x
         end where
      end if
   end function clip_values_1d
   !-------------------------------------------------------------------------------
   function clip_values_2d(x, min_x, max_x, varname, tolerance) result(error_message)
      real(r8), intent(inout) :: x(:,:)
      real(r8), intent(in) :: min_x
      real(r8), intent(in) :: max_x
      character(len=*), intent(in) :: varname
      real(r8), intent(in), optional :: tolerance
      real(r8) :: tolerance_local
      character(max_char_len) :: error_message

      error_message = ''
      tolerance_local = 0._r8
      if (present(tolerance)) then
         tolerance_local = tolerance
      end if

      ! look for values less than threshold
      if (any(x < min_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x < (min_x - tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x < (min_x - tolerance_local)), ' values below threshold ', &
                     '; min = ', minval(x)
         end if
         ! Clip values
         where (x < min_x)
            x = min_x
         endwhere
      end if

      ! Look for values greater than threshold
      if (any(x > max_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x > (max_x + tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x > (max_x + tolerance_local)), ' values above threshold ', &
                     '; max = ', maxval(x)
         end if
         ! Clip values
         where (x > max_x)
            x = max_x
         end where
      end if
   end function clip_values_2d
   !-------------------------------------------------------------------------------
   function clip_values_3d(x, min_x, max_x, varname, tolerance) result(error_message)
      real(r8), intent(inout) :: x(:,:,:)
      real(r8), intent(in) :: min_x
      real(r8), intent(in) :: max_x
      character(len=*), intent(in) :: varname
      real(r8), intent(in), optional :: tolerance
      real(r8) :: tolerance_local
      character(max_char_len) :: error_message

      error_message = ''
      tolerance_local = 0._r8
      if (present(tolerance)) then
         tolerance_local = tolerance
      end if

      ! look for values less than threshold
      if (any(x < min_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x < (min_x - tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x < (min_x - tolerance_local)), ' values below threshold ', &
                     '; min = ', minval(x)
         end if
         ! Clip values
         where (x < min_x)
            x = min_x
         endwhere
      end if

      ! Look for values greater than threshold
      if (any(x > max_x)) then
         ! Raise warning? Only if outside tolerance
         if (any(x > (max_x + tolerance_local))) then
            write(error_message,*) 'WARNING: ' // trim(varname) // ': ', &
                     count(x > (max_x + tolerance_local)), ' values above threshold ', &
                     '; max = ', maxval(x)
         end if
         ! Clip values
         where (x > max_x)
            x = max_x
         end where
      end if
   end function clip_values_3d
   !-------------------------------------------------------------------------------
   subroutine handle_error(error_message, fatal)
      use cam_abortutils, only: endrun
      character(len=*), intent(in) :: error_message
      logical, intent(in), optional :: fatal
      logical :: fatal_local = .true.

      ! Allow passing of an optional flag to not stop the run if an error is
      ! encountered. This allows this subroutine to be used when inquiring if a
      ! variable exists without failing.
      if (present(fatal)) then
         fatal_local = fatal
      else
         fatal_local = .true.
      end if

      ! If we encounter an error, fail if we require success. Otherwise do
      ! nothing and return silently.
      if (len(trim(error_message)) > 0) then
         if (fatal_local) then
            call endrun(trim(error_message))
         else
            print *, trim(error_message)
         end if
      end if
   end subroutine handle_error
   !----------------------------------------------------------------------------
end module radiation_utils
