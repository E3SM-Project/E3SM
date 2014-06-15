module interp

use control, only: verbose

implicit none
save
private

public :: interp_driver

contains

subroutine interp_driver(ncidi, ncido, f_in, f_out, vari, varo, nxi, &
                         nzi, nyi, nxo, nzo, nyo)

   ! first cut -- horizontal interpolation only
   ! Assume COARDS ordering, i.e., xyz

   use shr_kind_mod,     only: r8 => shr_kind_r8
   use fill_positions,   only: varspec_t
   use interpolate_data, only: lininterp, lininterp_init, interp_type, lininterp_finish

   implicit none

   include 'netcdf.inc'

   ! arguments

   integer, intent(in) :: ncidi, ncido
   integer, intent(in) :: nxi, nyi, nzi
   integer, intent(in) :: nxo, nyo, nzo

   type(varspec_t), intent(in) :: vari, varo

   real(r8), target, intent(in)  :: f_in(nxi,nyi,nzi)
   real(r8),         intent(out) :: f_out(nxo,nyo,nzo)

   ! Local workspace
   integer :: i,j,k
   integer :: numxis, numxin
   integer :: numxo
   integer :: jj, jjs, jjn
   integer :: count

   real(r8) :: x_in(nxi), x_out(nxo), y_in(nyi), z_in(nzi), z_out(nzo)
   real(r8), allocatable :: y_out(:)
   integer, parameter :: use_bnd_value_to_extrap = 1
   integer, parameter :: cyclic_bnds = 2
   type(interp_type) :: xwgts, ywgts, zwgts

   logical  :: do_vert_interp, do_horiz_interp
   real(r8), allocatable :: col_in(:), col_out(:)
   real(r8), pointer :: f_z(:,:,:)
   !-----------------------------------------------------------------------------

   ! Check whether vertical interpolation is needed.

   do_vert_interp = .false.

   if (nzi > 1) then

      ! If nzi > 1 then assume a vertical coordinate is present.  And if the input
      ! variable has a vertical coordinate then so does the output, even if its size
      ! is 1.
      call wrap_get_var_double(ncidi, vari%z_vid, z_in)
      call wrap_get_var_double(ncido, varo%z_vid, z_out)

      ! If number of input and output levels is the same, check whether the coordinates
      ! are different.  Otherwise no need to interpolate.
      if (nzi == nzo) then
         do k = 1, nzi
            if ( abs(z_in(k) - z_out(k)) > 1.e-4 ) do_vert_interp = .true.
         end do
      else
         ! If number of levels changes on output the interpolation is required.
         do_vert_interp = .true.
      end if

   end if

   if (do_vert_interp) then

      if (verbose) then
         write(6,*)'  doing vertical interpolation'
      end if

      allocate(f_z(nxi,nyi,nzo))
      allocate(col_in(nzi))
      allocate(col_out(nzo))

      call lininterp_init(z_in, nzi, z_out, nzo, use_bnd_value_to_extrap, zwgts)

      ! interpolate a column at a time
      do j = 1, nyi
         do i = 1, nxi
            do k = 1, nzi
               col_in(k) = f_in(i,j,k)
            end do
            call lininterp(col_in, nzi, col_out, nzo, zwgts)
            do k = 1, nzo
               f_z(i,j,k) = col_out(k)
            end do
         end do
      end do

      deallocate(col_in)
      deallocate(col_out)
   else
      f_z => f_in
   end if

   ! Get x and y coordinates
   call wrap_get_var_double(ncidi, vari%x_vid, x_in)
   call wrap_get_var_double(ncidi, vari%y_vid, y_in)
   call wrap_get_var_double(ncido, varo%x_vid, x_out)

   if (nyo > 1) then
      allocate(y_out(nyo))
   else
      ! This is the unstructured grid mode.  The number of global
      ! columns is nxo
      allocate(y_out(nxo))
   end if
   call wrap_get_var_double(ncido, varo%y_vid, y_out)

   ! Check whether horizontal interpolation is needed

   do_horiz_interp = .false.

   if (nxi == nxo  .and.  nyi == nyo) then

      ! If grids are the same shape then need to check coordinates
      do i = 1, nxi
         if ( abs(x_in(i) - x_out(i)) > 1.e-4 ) do_horiz_interp = .true.
      end do
      do j = 1, nyi
         if ( abs(y_in(j) - y_out(j)) > 1.e-4 ) do_horiz_interp = .true.
      end do

   else
      ! If grids aren't the same shape then interpolation is needed
      do_horiz_interp = .true.
   end if

   if (do_horiz_interp) then

      if (verbose) then
         write(6,*)'  doing horizontal interpolation'
      end if

      call lininterp_init(x_in, nxi, x_out, nxo, cyclic_bnds, xwgts)

      if (nyo > 1) then
         ! Output grid is rectangular.
         call lininterp_init(y_in, nyi, y_out, nyo, use_bnd_value_to_extrap, ywgts)
     
         do k = 1, nzo
            call lininterp(f_z(:,:,k), nxi, nyi, f_out(:,:,k), nxo, nyo, xwgts, ywgts)
         end do
      else	
         ! Output grid is unstructured.
         call lininterp_init(y_in, nyi, y_out, nxo, use_bnd_value_to_extrap, ywgts)

         do k = 1, nzo
            call lininterp(f_z(:,:,k), nxi, nyi, f_out(:,1,k), nxo, xwgts, ywgts)
         end do
      end if

      call lininterp_finish(xwgts)
      call lininterp_finish(ywgts)

   else
      f_out = f_z
   end if

   deallocate(y_out)
   if (do_vert_interp) then
      deallocate(f_z)
   end if

end subroutine interp_driver

end module interp
