program Test_Lib
	use pio
	use pio_kinds
	use piovdc
	implicit none
	
	include 'mpif.h'
	type (Var_desc_t)	:: var_handle, var_handle_no_comp
	character		:: vdf_path*100, vdc_path*100, var*40
	type (File_desc_t)	:: file_handle, file_handle_no_comp
	type (IOsystem_desc_t)	:: iosystem
	type (io_desc_t)	:: iodesc
	integer(i4)		:: rank, ierr, iostat,  dim_ids(3), nioprocs, nprocs, dims(3), dpp, n
	integer(i4),allocatable		:: compdof(:)
	integer(kind=PIO_OFFSET) :: iocount(3), iostart(3)
	real (r4),  allocatable :: array(:), read_array(:)
#ifdef 	DEBUG
	double precision	:: start, end
#endif


	vdf_path = '/ptmp/ypolius/Data/libbench.vdf' //CHAR(0)
	vdc_path = '/ptmp/ypolius/Data/benchdata.nc'
	var = 'vx' // CHAR(0)

#ifdef DEBUG 
	if(rank .eq. 0 ) then
	    print *, 'Calling piovdc_init...'
	endif
#endif
	nioprocs = 32
	call piovdc_init(vdf_path,  nioprocs, iostart, iocount, rank, ierr)


#ifdef DEBUG
	    print *, 'piovdc_init returned, rank: ' , rank, ' nioprocs: ' , nioprocs, ' iostart: ', iostart, ' iocount: ' , iocount
#endif


#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'Initiating PIO...'
	endif
#endif

	call PIO_init(rank, MPI_COMM_WORLD, nioprocs, nioprocs, 1, PIO_rearr_box, iosystem, 0)		

#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'PIO Initiated'
	endif
#endif
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	call piovdc_gen_vdc2_handles(file_handle, var_handle, iosystem, iocount,  vdf_path, ierr)

	dims = var_handle%dims

	dpp = product(dims) / nprocs !INT(temp_dpp)
	
	if(allocated(array)) then
		deallocate(array)
	endif

	allocate(array(dpp))
	allocate(read_array(dpp))

#ifdef DEBUG
	if(rank .eq. 0 ) then
		print *, 'Allocated write and read arrays'
		print *, 'File type', file_handle%iotype
	endif
#endif
	if(allocated(compdof)) then
	    deallocate(compdof)
	endif

	allocate(compdof(dpp))

#ifdef DEBUG	
	if(rank .eq. 0 ) then
	    print *, 'Allocated compdof'
	    print *, 'Filling compdof array...dpp-comp: ', dpp, ' dpp-io: ', product(dims)/nioprocs, ' dims: ' , dims, ' int limit: ' , huge(dims(1))
	endif
#endif

	do n = 1, dpp
	    compdof(n) =  n + rank * dpp !INT(REAL(n) + REAL(rank) * REAL(dpp))
	    array(n) = 53.0
		if(compdof(n) .lt. 0) then
			print *, ' n: ' , n , ' compdof(n): ' , compdof(n), ' array(n): ', array(n)
		end if
	end do

#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'Filled compdof array'
	    print *, 'Filled data array'
	endif
#endif
	
#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'Initializing decomposition...'
	endif
#endif

#ifdef DEBUG
	start = MPI_WTIME()
#endif
	if(PIO_iam_iotask(iosystem)) then

#ifdef DEBUG
print *, 'IAM IO TASK: ' , rank, ' iostart: ', iostart, ' iocount: ', iocount , ' dims: ', dims
#endif
		call PIO_initdecomp(iosystem, PIO_real, dims, compdof, iodesc, iostart, iocount)
	else
		call PIO_initdecomp(iosystem, PIO_real, dims, compdof, iodesc)
	endif

	print *, 'Decomp initialized'

#ifdef DEBUG
	end = MPI_WTIME()
	print *, 'Got time, writing darray'
#endif

#ifdef DEBUG
	if(rank .eq. 0) then 
	    print *, 'Decomposition initialized'
	    print *, 'Opening var for writing...(vx0)'
	endif
#endif

	call piovdc_open_var(0, var, -1, var_handle, ierr)

#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'Opened var for writing'
	endif
#endif

	call PIO_write_darray(file_handle, var_handle, iodesc,  array, iostat)

#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'Wrote data array'
	endif
#endif

#ifdef DEBUG
	print *, 'Rank: ', rank, ' vdc write time: ', MPI_WTIME() - end
	print *, 'Rank: ', rank, 'Decomposition rearrangment runtime: ' , end - start
#endif

#ifdef DEBUG
	print *, 'Rank: ', rank, 'Aggregate runtime: ' , MPI_WTIME() - start
#endif

#ifdef DEBUG
	if(rank .eq. 0) then
		print *, 'Attempting to read back data'
	endif
#endif
#ifdef DEBUG
	start = MPI_WTIME()
#endif
	call PIO_read_darray(file_handle, var_handle, iodesc,  read_array, iostat)

#ifdef DEBUG
	print *, 'Rank: ', rank, ' vdc read time: ' , MPI_WTIME() - start
#endif
	ierr = PIO_CreateFile(iosystem, file_handle_no_comp, PIO_iotype_pnetcdf, vdc_path, PIO_clobber)
	var_handle_no_comp%type = PIO_iotype_pnetcdf
	iostat = PIO_def_dim(file_handle_no_comp, 'z', dims(3), dim_ids(3))
	iostat = PIO_def_dim(file_handle_no_comp, 'y', dims(2), dim_ids(2))
	iostat = PIO_def_dim(file_handle_no_comp, 'x', dims(1), dim_ids(1))
	iostat = PIO_def_var(file_handle_no_comp, 'vx', PIO_real, dim_ids, var_handle_no_comp)
	ierr = PIO_enddef(file_handle_no_comp)

#ifdef DEBUG
	start = MPI_WTIME()
#endif

	call pio_write_darray(file_handle_no_comp, var_handle_no_comp, iodesc, array, iostat)

#ifdef DEBUG
	if(rank .eq. 0) then
		print *, 'Pure NC write_darray runtime: ' , MPI_WTIME() - start
	endif
#endif

#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'Wrote data array'
	end if
#endif

	!clean up PIO and MPI
	call PIO_CloseFile(file_handle_no_comp)

#ifdef DEBUG
	if(rank .eq. 0 ) then
		print *, 'Closed PIO file'
	endif
#endif

	call PIO_Finalize(iosystem, ierr)

#ifdef DEBUG
	if(rank .eq. 0 ) then
		print *, 'Finalized PIO'	
	endif
#endif

	deallocate(compdof)
	deallocate(array)

	!force non-io tasks to barrier until io is complete

	call MPI_Finalize(ierr)
	stop
endprogram
