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
      module procedure clip_values_1d, clip_values_2d
   end interface clip_values

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

   subroutine calculate_heating_rate(fluxes, pint, heating_rate)

      use physconst, only: gravit
      use mo_fluxes_byband, only: ty_fluxes_byband

      ! Inputs
      type(ty_fluxes_byband), intent(in) :: fluxes
      real(r8), intent(in) :: pint(:,:)

      ! Output heating rate; same size as pint with one fewer vertical level
      real(r8), intent(out) :: heating_rate(:,:)

      ! Loop indices
      integer :: icol, ilev

      ! Everyone needs a name
      character(len=32) :: subname = 'calculate_heating_rate'

      ! Get dimension sizes and make sure arrays conform
      call assert(size(pint,1) == size(fluxes%flux_up,1), subname // ': sizes do not conform.')
      call assert(size(pint,2) == size(fluxes%flux_up,2), subname // ': sizes do not conform.')
      call assert(size(heating_rate,1) == size(fluxes%flux_up,1), subname // ': sizes do not conform.')
      call assert(size(heating_rate,2) == size(fluxes%flux_up,2)-1, subname // ': sizes do not conform.')

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
               fluxes%flux_up(icol,ilev+1) - fluxes%flux_up(icol,ilev) - &
               fluxes%flux_dn(icol,ilev+1) + fluxes%flux_dn(icol,ilev) &
            ) * gravit / (pint(icol,ilev+1) - pint(icol,ilev))
         end do
      end do

   end subroutine calculate_heating_rate

   !----------------------------------------------------------------------------

   subroutine clip_values_1d(x, min_x, max_x, varname, warn)
      real(r8), intent(inout) :: x(:)
      real(r8), intent(in) :: min_x
      real(r8), intent(in) :: max_x
      character(len=*), intent(in), optional :: varname
      logical, intent(in), optional :: warn

      logical :: warn_local

      warn_local = .false.
      if (present(warn)) then
         warn_local = warn
      end if

      ! Look for values less than threshold
      if (any(x < min_x)) then
         ! Raise warning?
         if (warn_local) then
            if (present(varname)) then
               print *, module_name // ' warning: ', &
                        count(x < min_x), ' values are below threshold for variable ', &
                        trim(varname), '; min = ', minval(x)
            else
               print *, module_name // ' warning: ', &
                        count(x < min_x), ' values are below threshold; min = ', minval(x)
            end if
         end if

         ! Clip values
         where (x < min_x)
            x = min_x
         endwhere
      end if

      ! Look for values greater than threshold 
      if (any(x > max_x)) then 
         ! Raise warning?
         if (warn_local) then
            if (present(varname)) then
               print *, module_name // ' warning: ', &
                        count(x > max_x), ' values are above threshold for variable ', &
                        trim(varname), '; max = ', maxval(x)
            else
               print *, module_name // ' warning: ', &
                        count(x > max_x), ' values are above threshold; max = ', maxval(x)
            end if
         end if

         ! Clip values
         where (x > max_x)
            x = max_x
         end where
      end if
   end subroutine clip_values_1d

   subroutine clip_values_2d(x, min_x, max_x, varname, warn)
      real(r8), intent(inout) :: x(:,:)
      real(r8), intent(in) :: min_x
      real(r8), intent(in) :: max_x
      character(len=*), intent(in), optional :: varname
      logical, intent(in), optional :: warn

      logical :: warn_local

      warn_local = .false.
      if (present(warn)) then
         warn_local = warn
      end if

      ! look for values less than threshold
      if (any(x < min_x)) then
         ! Raise warning?
         if (warn_local) then
            if (present(varname)) then
               print *, module_name // ' warning: ', &
                        count(x < min_x), ' values are below threshold for variable ', &
                        trim(varname), '; min = ', minval(x)
            else
               print *, module_name // ' warning: ', &
                        count(x < min_x), ' values are below threshold; min = ', minval(x)
            end if
         end if

         ! Clip values
         where (x < min_x)
            x = min_x
         endwhere
      end if

      ! Look for values greater than threshold
      if (any(x > max_x)) then 
         ! Raise warning?
         if (warn_local) then
            if (present(varname)) then
               print *, module_name // ' warning: ', &
                        count(x > max_x), ' values are above threshold for variable ', &
                        trim(varname), '; max = ', maxval(x)
            else
               print *, module_name // ' warning: ', &
                        count(x > max_x), ' values are above threshold; max = ', maxval(x)
            end if
         end if

         ! Clip values
         where (x > max_x)
            x = max_x
         end where
      end if

   end subroutine clip_values_2d

   !----------------------------------------------------------------------------

   subroutine handle_error(error_message, stop_on_error)
      use cam_abortutils, only: endrun
      character(len=*), intent(in) :: error_message
      logical, intent(in), optional :: stop_on_error
      logical :: stop_on_error_local = .true.

      ! Allow passing of an optional flag to not stop the run if an error is
      ! encountered. This allows this subroutine to be used when inquiring if a
      ! variable exists without failing.
      if (present(stop_on_error)) then
         stop_on_error_local = stop_on_error
      else
         stop_on_error_local = .true.
      end if

      ! If we encounter an error, fail if we require success. Otherwise do
      ! nothing and return silently.
      if (len(trim(error_message)) > 0) then
         if (stop_on_error_local) then
            call endrun(module_name // ': ' // error_message)
         end if
      end if
   end subroutine handle_error

   !----------------------------------------------------------------------------

end module radiation_utils
