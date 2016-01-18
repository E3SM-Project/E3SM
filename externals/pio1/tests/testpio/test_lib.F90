program Test_Lib
	use pio !everything is forwarded from the PIO module, should be the only use necessary
	use pio_kinds
	implicit none
	
	include 'mpif.h'
	type (Var_desc_t)	:: var_handle_no_comp !type handle for normal, uncompressed PIO variables
	type (Var_desc_t)	:: var_handle !type handle for compressed, VDC variables
	character		:: vdf_path*100 !location to save the vdf file and it's related data
	character		:: binary_path*100 !location to save the binary data
	type (File_desc_t)	:: file_handle !each open file requires a separate file handle, this one is for VDC data
	type (File_desc_t)	:: file_handle_no_comp !file handle for normal PIO variables
	type (IOsystem_desc_t)	:: iosystem !PIO type handle to hold PIO-specific information about a file's IO settings
	type (io_desc_t)	:: iodesc !PIO type handle to hold PIO-specific information about a file's IO decomposition
	integer(i4)		:: rank !MPI Rank of the current process
	integer(i4)		:: ierr !general error code variable
	integer(i4)		:: iostat !PIO-specific error code variable
	integer(i4)		:: dim_ids(3) !used in the uncompressed PIO for defining dimensions used in a *cdf file
	integer(i4)		:: nioprocs !used to tell PIO how many IO procs you want to use, functions as the max # of IO procs wanted when using compression, less may be used
	integer(i4)		:: nprocs !the # of processes involve in the MPI run
	integer(i4)		:: dims(3) !the 3D grid size used to write VDC data
	integer(i4)		:: n !counter
	integer(kind=PIO_Offset) :: dpp !data per process, the amount of data each MPI task contributes to the overall file
	integer(kind=PIO_Offset),allocatable :: compdof(:) !computational degrees-of-freedom, this array holds the mapping from the local 
														!slice of computational data to the global grid
	real (r4),  allocatable :: array(:), read_array(:) !arrays holding the local computational data
#ifdef 	DEBUG
	double precision	:: start, end, temp !timing variables
#endif
	
	!first set locals for vdc compression and the uncompressed data path
	dims = (/1024, 1024, 1024/)
	vdf_path = '/glade/scratch/ypolius/piovdc/libbench.vdf'
	binary_path = '/glade/scratch/ypolius/piovdc/benchdata.nc'
	nioprocs = 64
	
	!init MPI and retrieve MPI-specific info
	call MPI_init(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'Initiating PIO...'
	endif
#endif
	
	!call PIO_init to initiate iosystem
	call PIO_init(rank, MPI_COMM_WORLD, nioprocs, nioprocs, 1, PIO_rearr_box, iosystem, 0)		
	
#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'PIO Initiated procs: ', nioprocs
	endif
#endif

	!set data-per-process to be the # of grid elements / # of computational procs
	!conversions to PIO_Offset int to allow for extremely large dims, ex 2048 x 2048 x 2048
	dpp = int(dims(1), kind=PIO_Offset) * int(dims(2), kind=PIO_Offset) * int(dims(3), kind=PIO_Offset) / int(nprocs, kind=PIO_Offset)
	
	!allocate local memory arrays	
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

	!allocate compdof
	if(allocated(compdof)) then
	    deallocate(compdof)
	endif

	allocate(compdof(dpp))

#ifdef DEBUG	
	if(rank .eq. 0 ) then
	    print *, 'Allocated compdof'
	    print *, 'Filling compdof array...dpp-comp: ', dpp, ' dpp-io: ', product(dims)/nioprocs, ' dims: ' , dims, ' int limit: ' , huge(compdof(1)), ' sample calc: ' , product(dims)
	endif
#endif

	!setup mapping from local data to global grid, fill local data with sample data
	!sample data (array) is a constant for this example, but can be any floating point value so user data may be used for this
	!this example uses a simple linear mapping, ex rank 0 points to global data [0, dpp], rank 1 points to global data[dpp, dpp *2]
	!any mapping is possible, care must be taken to make sure that the mapping is strictly 1-1 or PIO will report errors
	do n = 1, dpp
	    compdof(n) =  int(n + rank * dpp,kind=pio_offset) !INT(REAL(n) + REAL(rank) * REAL(dpp))
	    array(n) = 53.0
		if(compdof(n) .lt. 0) then
			print *, ' n: ' , n , ' compdof(n): ' , compdof(n), ' array(n): ', array(n)
		end if
	end do

#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'Filled compdof array: ', compdof(1) , '-', compdof(2)
	    print *, 'Filled data array: '!, array(1), '-' , array(2)
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

	!call init_decomp in order to setup the IO decomposition with PIO
        ! The optional parameter num_ts is required to indicate the vdc output method
	call PIO_initdecomp(iosystem, PIO_real, dims, compdof, iodesc, num_ts=10) 
	
	!example using optional bsize and # timesteps specifiers
	!call PIO_initdecomp(iosystem, PIO_real, dims, compdof, iodesc, bsize=(/128, 128, 128/), num_ts=30)		
#ifdef DEBUG
	if(rank .eq. 0) then 
	    print *, 'Decomposition initialized'
	    print *, 'Creating vdf file'
	endif
        print *, 'Rank: ', rank, 'Decomposition rearrangment runtime: ' , MPI_WTIME() - start
#endif

	!use create file with VDC2 io type to begin the creation of a VDF metadata file, not valid until enddef is called
	ierr = PIO_CreateFile(iosystem, file_handle, PIO_iotype_vdc2, vdf_path, PIO_clobber)

#ifdef DEBUG
	if(rank .eq. 0) then 
	    print *, 'VDF file created'
	    print *, 'Opening vdf var for writing...(vx0)'
	endif
#endif
	
	!define the variables that will be written into the VDC 
	!VDC WRITING DOES NOT REQUIRE CREATING DIMS, THERE ARE ALWAYS 3, dims = (/x, y, z/)
	iostat = PIO_def_var(file_handle, 'vx' , PIO_real,var_handle)

#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'Opened var for writing'
	endif
#endif

	!finally call enddef to have the VDF metadata file written out
	ierr = PIO_enddef(file_handle)
	call PIO_setframe(var_handle, int(0, PIO_OFFSET))
#ifdef DEBUG
	if(rank .eq. 0 ) then
	    print *, 'Ended VDF definition'
	endif
	start = MPI_WTIME()
#endif

	!to write data call PIO_write_darray, the only difference with compressed vs uncompressed
	!writing is that compressed writing requires that the the user inputs the current time step 
	!corresponding to the variable about to be written
	call PIO_write_darray(file_handle, var_handle, iodesc,  array, iostat)

#ifdef DEBUG
	print *, 'Rank: ', rank, ' vdc write time: ', MPI_WTIME() - start
#endif

#ifdef DEBUG
	if(rank .eq. 0) then
		print *, 'Attempting to read back data'
	endif
	start = MPI_WTIME()
#endif


	!to read data call PIO_read_darray, the only difference with compressed vs uncompressed
	!reading is that compressed reading requires that the the user inputs the current time step 
	!corresponding to the variable about to be read
	call PIO_read_darray(file_handle, var_handle, iodesc,  read_array, iostat)

#ifdef DEBUG
	print *, 'Rank: ', rank, ' vdc read time: ' , MPI_WTIME() - start
#endif

	!Setup for UNCOMPRESSED files
	
	!Same call as the VDF file, but we switch the IO type to netcdf
	ierr = PIO_CreateFile(iosystem, file_handle_no_comp, PIO_iotype_pnetcdf, binary_path, PIO_clobber)
	
	!define the dimensions to be used with the file PIO writes
	iostat = PIO_def_dim(file_handle_no_comp, 'z', dims(3), dim_ids(3))
	iostat = PIO_def_dim(file_handle_no_comp, 'y', dims(2), dim_ids(2))
	iostat = PIO_def_dim(file_handle_no_comp, 'x', dims(1), dim_ids(1))
	
	!define variables to be written
	iostat = PIO_def_var(file_handle_no_comp, 'vx', PIO_real, dim_ids, var_handle_no_comp)
	
	!end definition to make file valid
	ierr = PIO_enddef(file_handle_no_comp)

#ifdef DEBUG
	start = MPI_WTIME()
#endif

	!write data, like compressed files but no timestep required at the end
	call pio_write_darray(file_handle_no_comp, var_handle_no_comp, iodesc, array, iostat)

#ifdef DEBUG

	print *, 'Rank: ', rank , ' pure NC write_darray runtime: ' , MPI_WTIME() - start

#endif

	!close non-compressed files, compressed files are automatically closed
	call PIO_CloseFile(file_handle_no_comp)

#ifdef DEBUG
	if(rank .eq. 0 ) then
		print *, 'Closed PIO file'
	endif
#endif

	!clean up PIO
	call PIO_Finalize(iosystem, ierr)

#ifdef DEBUG
	if(rank .eq. 0 ) then
		print *, 'Finalized PIO'	
	endif
#endif

	!clean up data arrays
	deallocate(compdof)
	deallocate(array)

	!clean up MPI

	call MPI_Finalize(ierr)
	stop
endprogram
