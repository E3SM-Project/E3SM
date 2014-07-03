
      module mo_inter

      use shr_kind_mod, only : r8 => shr_kind_r8
      use cam_logfile,  only : iulog

      implicit none

      integer :: nintervals
      integer, allocatable :: xi(:), xcnt(:)
      real(r8), allocatable    :: xfrac(:,:)

      contains

      subroutine inter2( ng, xg, yg, n, x, y, ierr )
!-----------------------------------------------------------------------------
!   purpose:
!   map input data given on single, discrete points onto a set of target
!   bins.
!   the original input data are given on single, discrete points of an
!   arbitrary grid and are being linearly interpolated onto a specified set
!   of target bins.  in general, this is the case for most of the weighting
!   functions (action spectra, molecular cross section, and quantum yield
!   data), which have to be matched onto the specified wavelength intervals.
!   the average value in each target bin is found by averaging the trapezoi-
!   dal area underneath the input data curve (constructed by linearly connec-
!   ting the discrete input values).
!   some caution should be used near the endpoints of the grids.  if the
!   input data set does not span the range of the target grid, an error
!   message is printed and the execution is stopped, as extrapolation of the
!   data is not permitted.
!   if the input data does not encompass the target grid, use addpnt to
!   expand the input array.
!-----------------------------------------------------------------------------
!   parameters:
!   ng  - integer, number of bins + 1 in the target grid                  (i)
!   xg  - real, target grid (e.g., wavelength grid);  bin i is defined    (i)
!         as [xg(i),xg(i+1)] (i = 1..ng-1)
!   yg  - real, y-data re-gridded onto xg, yg(i) specifies the value for  (o)
!         bin i (i = 1..ng-1)
!   n   - integer, number of points in input grid                         (i)
!   x   - real, grid on which input data are defined                      (i)
!   y   - real, input y-data                                              (i)
!-----------------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)  :: ng, n
      integer, intent(out) :: ierr
      real(r8), intent(in)     :: x(n)
      real(r8), intent(in)     :: y(n)
      real(r8), intent(in)     :: xg(ng)
      real(r8), intent(out)    :: yg(ng)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer :: ngintv
      integer :: i, k, jstart
      real(r8)    :: area, xgl, xgu
      real(r8)    :: darea, slope
      real(r8)    :: a1, a2, b1, b2

      ierr = 0
!-----------------------------------------------------------------------------
!  	... test for correct ordering of data, by increasing value of x
!-----------------------------------------------------------------------------
      do i = 2,n
         if( x(i) <= x(i-1) ) then
            ierr = 1
            write(iulog,*) 'inter2: x coord not monotonically increasing'
            return
         end if
      end do

      do i = 2,ng
        if( xg(i) <= xg(i-1) ) then
           ierr = 2
           write(iulog,*) 'inter2: xg coord not monotonically increasing'
           return
        end if
      end do

!-----------------------------------------------------------------------------
! 	... check for xg-values outside the x-range
!-----------------------------------------------------------------------------
      if( x(1) > xg(1) .or. x(n) < xg(ng) ) then
          write(iulog,*) 'inter2: data does not span grid'
          write(iulog,*) '        use addpnt to expand data and re-run'
!         call endrun
      end if

!-----------------------------------------------------------------------------
!  	... find the integral of each grid interval and use this to 
!           calculate the average y value for the interval      
!           xgl and xgu are the lower and upper limits of the grid interval
!-----------------------------------------------------------------------------
      jstart = 1
      ngintv = ng - 1
      do i = 1,ngintv
!-----------------------------------------------------------------------------
! 	... initalize:
!-----------------------------------------------------------------------------
         area = 0._r8
         xgl = xg(i)
         xgu = xg(i+1)
!-----------------------------------------------------------------------------
!  	... discard data before the first grid interval and after the 
!           last grid interval
!           for internal grid intervals, start calculating area by interpolating
!           between the last point which lies in the previous interval and the
!           first point inside the current interval
!-----------------------------------------------------------------------------
         k = jstart
         if( k <= n-1 ) then
!-----------------------------------------------------------------------------
!  	... if both points are before the first grid, go to the next point
!-----------------------------------------------------------------------------
            do
               if( x(k+1)  <= xgl ) then
                  jstart = k - 1
                  k = k+1
                  if( k <= n-1 ) then
		     cycle
		  else
		     exit
		  end if
	       else
		  exit
               end if
	    end do
!-----------------------------------------------------------------------------
!  	... if the last point is beyond the end of the grid,
!           complete and go to the next grid
!-----------------------------------------------------------------------------
	    do
               if( k <= n-1 .and. x(k) < xgu ) then          
                  jstart = k-1
!-----------------------------------------------------------------------------
! 	... compute x-coordinates of increment
!-----------------------------------------------------------------------------
                  a1 = max( x(k),xgl )
                  a2 = min( x(k+1),xgu )
!-----------------------------------------------------------------------------
!  	... if points coincide, contribution is zero
!-----------------------------------------------------------------------------
                  if( x(k+1) == x(k) ) then
                     darea = 0._r8
                  else
                     slope = (y(k+1) - y(k))/(x(k+1) - x(k))
                     b1    = y(k) + slope*(a1 - x(k))
                     b2    = y(k) + slope*(a2 - x(k))
                     darea = .5_r8*(a2 - a1)*(b2 + b1)
                  end if
!-----------------------------------------------------------------------------
!  	... find the area under the trapezoid from a1 to a2
!-----------------------------------------------------------------------------
                  area = area + darea
                  k = k+1
                  cycle
	       else
		  exit
               end if
	    end do
         end if
!-----------------------------------------------------------------------------
!  	... calculate the average y after summing the areas in the interval
!-----------------------------------------------------------------------------
         yg(i) = area/(xgu - xgl)
      end do

      end subroutine inter2

      subroutine inter_inti( ng, xg, n, x )
!-----------------------------------------------------------------------------
! 	... initialization
!-----------------------------------------------------------------------------

      use abortutils,   only : endrun

      implicit none
      
!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in) :: ng, n
      real(r8), intent(in)    :: xg(ng)
      real(r8), intent(in)    :: x(n)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer :: i, ii, iil, astat
      integer :: ndim(1)
#ifdef DEBUG
      write(iulog,*) 'inter3: diagnostics; ng,n = ',ng,n
      write(iulog,'('' xg '' )')
      write(iulog,'(1p,5e21.13)') xg(:ng)
      write(iulog,'('' x  '' )')
      write(iulog,'(1p,5e21.13)') x(:n)
#endif
      allocate( xi(ng), xcnt(ng-1), stat=astat )
      if( astat /= 0 ) then
	 write(iulog,*) 'inter_inti: failed to allocate wrk arrays; error = ',astat
	 call endrun
      else
	 xi(:)   = 0
	 xcnt(:) = 0
      end if
      iil = 1
      do i = 1,ng
         do ii = iil,n-1
            if( xg(i) < x(ii) ) then
               xi(i) = ii - 1
	       iil   = ii
	       exit
            end if
         end do
      end do
      nintervals = count( xi(:) /= 0 )
      if( nintervals == 0 ) then
         write(iulog,*) 'inter_inti: wavelength grids do not overlap'
	 call endrun
      else
	 nintervals = nintervals - 1
      end if
      xcnt(1:nintervals) = xi(2:nintervals+1) - xi(1:nintervals) + 1
      ndim(:) = maxval( xcnt(1:nintervals) )
      allocate( xfrac(ndim(1),nintervals),stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'inter_inti: failed to allocate wrk array; error = ',astat
	 call endrun
      else
         xfrac(:,:) = 1._r8
      end if
      do i = 1,nintervals
        iil = xi(i)
        xfrac(1,i) = (min( x(iil+1),xg(i+1) ) - xg(i))/(x(iil+1) - x(iil))
        if( xcnt(i) > 1 ) then
           iil = xi(i) + xcnt(i) - 1
           xfrac(xcnt(i),i) = (xg(i+1) - x(iil))/(x(iil+1) - x(iil))
        end if
      end do
      write(iulog,*) 'inter_inti: diagnostics; ng,n,nintervals = ',ng,n,nintervals
      write(iulog,'('' xi'')')
      write(iulog,'(10i4)') xi(1:nintervals+1)
      write(iulog,'('' xcnt'')')
      write(iulog,'(10i4)') xcnt(1:nintervals)
      write(iulog,'('' xfrac'')')
      do i = 1,nintervals
         write(iulog,'(1p,5e21.13)') xfrac(1:xcnt(i),i)
      end do

      end subroutine inter_inti

      subroutine inter3( ng, xg, yg, n, x, y )
!-----------------------------------------------------------------------------
!   purpose:
!   map input data given on a set of bins onto a different set of target
!   bins.
!   the input data are given on a set of bins (representing the integral
!   of the input quantity over the range of each bin) and are being matched
!   onto another set of bins (target grid).  a typical example would be an
!   input data set spcifying the extra-terrestrial flux on wavelength inter-
!   vals, that has to be matched onto the working wavelength grid.
!   the resulting area in a given bin of the target grid is calculated by
!   simply adding all fractional areas of the input data that cover that
!   particular target bin.
!   some caution should be used near the endpoints of the grids.  if the
!   input data do not span the full range of the target grid, the area in
!   the "missing" bins will be assumed to be zero.  if the input data extend
!   beyond the upper limit of the target grid, the user has the option to
!   integrate the "overhang" data and fold the remaining area back into the
!   last target bin.  using this option is recommended when re-gridding
!   vertical profiles that directly affect the total optical depth of the
!   model atmosphere.
!-----------------------------------------------------------------------------
!   parameters:
!   ng     - integer, number of bins + 1 in the target grid               (i)
!   xg     - real, target grid (e.g. working wavelength grid);  bin i     (i)
!            is defined as [xg(i),xg(i+1)] (i = 1..ng-1)
!   yg     - real, y-data re-gridded onto xg;  yg(i) specifies the        (o)
!            y-value for bin i (i = 1..ng-1)
!   n      - integer, number of bins + 1 in the input grid                (i)
!   x      - real, input grid (e.g. data wavelength grid);  bin i is      (i)
!            defined as [x(i),x(i+1)] (i = 1..n-1)
!   y      - real, input y-data on grid x;  y(i) specifies the            (i)
!            y-value for bin i (i = 1..n-1)
!-----------------------------------------------------------------------------

      implicit none
      
!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in) :: n, ng
      real(r8), intent(in)    :: xg(ng)
      real(r8), intent(in)    :: x(n)
      real(r8), intent(in)    :: y(n)
      real(r8), intent(out)   :: yg(ng)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer :: i, ii, iil

!-----------------------------------------------------------------------------
! 	... do interpolation
!-----------------------------------------------------------------------------
      yg(:) = 0._r8
      do i = 1,nintervals
	 iil = xi(i)
	 ii = xcnt(i)
	 if( ii == 1 ) then
	    yg(i) = xfrac(1,i)*y(iil)
	 else
	    yg(i) = dot_product( xfrac(1:ii,i),y(iil:iil+ii-1) )
	 end if
      end do
#ifdef DEBUG
      write(iulog,'('' y '')')
      write(iulog,'(1p,5e21.13)') y
      write(iulog,'('' yg '')')
      write(iulog,'(1p,5e21.13)') yg(1:nintervals)
#endif

      end subroutine inter3

      end module mo_inter
