!> @file libpiovdc.F90
!> @author  Yannick Polius <ypolius@ucar.edu>
!> @version 1.0
!> @date 03/22/2012
!> @brief The piovdc library for writing Vapor Data Collection (VDC) 2 data files
!> 
!> <br>
!> @details The piovdc library is used to write VDC2 data files in a	 
!> parallel manner using PIO. After the prerequisite library functions are 
!> used, a call to pio_writedarray is made, writing the passed 
!> data to an on disk VDC2 collection.<br>
!> PRE-REQUISITES: <br>
!> 	VDF meta-file must be generated, using either rawtovdf or vdfcreate
!>	if advanced features (wavelet type, compression ratios, or boundary type) are needed
!> 	VDF file requires VDC version to be 2, and requires the Waveletname,
!>	WaveletBoundaryMode, CompressionRations, and NumTransforms to be set.<br>
!> POST-EFFECTS: <br>
!>	After a successful write, VDC2 data will be in a directory located in
!>	the same directory as the vdf file, using the vdf name, appended with _data
!> 	(ex. ghost.vdf generates VDC2 data in the dir ghost_data in the vdf dir)
!>	If no compression is enabled, a single, uncompressed file will be 
!>	generated using PIO instead of a VDC
module piovdc
	use pio_kinds, only : i4, r4, pio_offset
	implicit none
	integer (i4)	:: vdc_dims(3), vdc_bsize(3), vdc_ts
	integer (kind=PIO_OFFSET)  :: vdc_iostart(3), vdc_iocount(3)	
contains

!> @brief subroutine checks start/count for out of bounds, adjusts if the start/count is too high, zeroes start if it is invalid
!> POST-EFFECTS:
!>	<br>all start/counts are now legal, non-IO tasks have zeroed start counts
!> @param[in] global_dims int(3) global grid dimensions
!> @param[in] rank int rank of current MPI task
!> @param[inout] start int(3) current MPI task global start
!> @param[inout] count int(3) current MPI task global count
subroutine adjust_bounds(global_dims, start, count, rank)
	real (r4), dimension(:), intent(in) :: global_dims
	integer(i4), intent(in) :: rank
	integer (kind=PIO_OFFSET), dimension(:), intent(inout) :: start, count

	!first check to ensure the start is legal

	if (start(1) .GT. global_dims(1) .OR. start(2) .GT. global_dims(2) .OR. start(3) .GT. global_dims(3)) then !outside of global bounds!

	   !negate everything, they're useless
#ifdef DEBUG
	   print *, ' rank: ' , rank, ' start: ' , start, ' count: ' , count , ' negated'
#endif
	   start = (/ 0, 0, 0/)
	   count = (/ 0, 0, 0/)
	else 
	   !start is legit but count might not be, check & adjust to the boundaries
	   if(count(1) + start(1) - 1 .GT. global_dims(1)) then
	      count(1) = global_dims(1) - start(1) + 1

	   endif
	   if(count(2) + start(2) - 1 .GT. global_dims(2)) then
	      count(2) = global_dims(2) - start(2) + 1

	   endif
	   if(count(3) + start(3) - 1 .GT. global_dims(3)) then
	      count(3) = global_dims(3) - start(3) + 1

	   endif
	end if
end subroutine

!> @brief subroutine that, given a global grid, VDC blocksize, and max # of nioprocs, will
!> automatically create an VDC optimized IO decomposition that uses the most possible IO tasks
!
!> POST-EFFECTS:
!>	<br>Each MPI Task is now either and IO task or a computational task. IO tasks have nonzero start/counts
!> @param[in] rank int rank of the current MPI task
!> @param[inout] nioprocs int represents the max possible # of IO procs, algorithm will try to get as close as possible to this # and return it in nioprocs
!> @param[in] blockdims int(3) global grid dimensions represented as VDC blocks
!> @param[out] start int(3) iostart for the current MPI task
!> @param[out] count int(3) iocount for the current MPI task
!> @param[in] bsize int(3) VDC block size
subroutine auto_get_start_count(rank, nioprocs, block_dims, start, count, bsize)
  use pio_kinds
  integer (kind=PIO_OFFSET), intent(out):: start(3), count(3)
  integer(i4), dimension(:), intent(in) :: bsize
  integer (i4), intent(in) :: rank
  real (r4), dimension(:), intent(in) 	:: block_dims

  integer (i4), intent(inout) :: nioprocs
  !locals
  real (r4)             :: proc_count
  integer (i4)             :: lpp, spp0, spp1, counter, slab_counter, calc_procs, nslabs, nlinesPslab
  logical		:: found

  found = .FALSE.
  nlinesPslab = CEILING(block_dims(2)) !max # of possible lines per slab PER TASK
  nslabs = CEILING(block_dims(3))
  calc_procs = -1
  if (nioprocs .EQ. 1) then
	nioprocs = 1
	start = (/0, 0, 0/)
	count = block_dims * bsize
  else
	do slab_counter=1, nslabs
	  do counter=1, nlinesPslab 
		proc_count =  CEILING(nlinesPslab / REAL(counter)) * CEILING(nslabs / REAL(slab_counter))
		!test to see if counter # of lines per processor per slab is possible
		if (nioprocs >= proc_count) then
			if (proc_count .gt. calc_procs) then
				calc_procs = proc_count ! return the actual # of io procs used
				count = (/ INT(block_dims(1) * bsize(1)), counter * bsize(2), slab_counter *bsize(3) /)
				start = (/ 0, mod(rank, INT(CEILING(nlinesPslab / REAL(counter)))) * counter * bsize(2), INT(rank / CEILING(nlinesPslab / REAL(counter))) * slab_counter * bsize(3)/) + 1
				call adjust_bounds(block_dims * bsize, start, count, rank)

				if(proc_count .eq. nioprocs) then !using max #of procs, suitable solution found (for now)
					found = .TRUE.
					exit				
				end if
				  if (found) then
					exit
				  end if
			end if
		end if
#ifdef DEBUG

#endif
	  end do
	  if (found) then
		exit
	  end if
	end do
  end if
 nioprocs = calc_procs
end subroutine

!> @brief subroutine that prepares the global grid to be split by the auto_start_count routine
!
!> POST-EFFECTS:
!>	<br> A valid IO decomposition is created that can be used with PIO
!> @param[in] rank int rank of the current
!> @param[in] data_dims int(3) size of the global grid
!> @param[in] vdc_bsize int(3) VDC block size
!> @param[out] iostart int(3) IO start for the current MPI task
!> @param[out] iocount int(3) IO count for the current MPI task
!> @param[inout] ioprocs int max # of IO procs, gets returned as the actual # used
subroutine init_vdc2(rank, data_dims, vdc_bsize, iostart, iocount, ioprocs)
  integer (kind=PIO_OFFSET), intent(out)  :: iostart(3), iocount(3)
  integer (i4), intent(in) :: rank
  integer(i4), dimension(:), intent(in) :: data_dims, vdc_bsize
  integer (i4), intent(inout):: ioprocs
  !locals
  real(r4) :: vdc_blocks(3)   
  integer (i4)	:: ierr

!#ifdef DEBUG
!  if(rank .eq. 0) then
	print *, 'Calling get start count...block_dims: ', data_dims/real(vdc_bsize), ' bsize: ' , vdc_bsize, ' ioprocs: ', ioprocs, ' dims: ', data_dims, ' rank: ',rank
!  endif
!#endif
  vdc_blocks = data_dims/real(vdc_bsize)
  
  call auto_get_start_count (rank, ioprocs, vdc_blocks, iostart, iocount, vdc_bsize)
  
!#ifdef DEBUG 
		print *, 'Retrieved VDF start count', iostart, '-', iocount, 'rank: ' , rank
!#endif

endsubroutine

end module piovdc
