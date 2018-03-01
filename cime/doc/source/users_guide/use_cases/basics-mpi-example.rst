.. _mpi-example:

MPI example
-----------

It is usually very helpful to assure that you can run a basic mpi parallel program on your machine prior to attempting a CIME port. 
Understanding how to compile and run the program fhello_world_mpi.F90 shown here could potentially save many hours of frustration.
::

   program fhello_world_mpi.F90
     use mpi
     implicit none
     integer ( kind = 4 ) error
     integer ( kind = 4 ) id
     integer p
     character(len=MPI_MAX_PROCESSOR_NAME) :: name
     integer clen
     integer, allocatable :: mype(:)
     real ( kind = 8 ) wtime

     call MPI_Init ( error )
     call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
     call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
     if ( id == 0 ) then
        wtime = MPI_Wtime ( )
     
	write ( *, '(a)' ) ' '
	write ( *, '(a)' ) 'HELLO_MPI - Master process:'
	write ( *, '(a)' ) '  FORTRAN90/MPI version'
	write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  An MPI test program.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  The number of processes is ', p
        write ( *, '(a)' ) ' '
     end if
     call MPI_GET_PROCESSOR_NAME(NAME, CLEN, ERROR)
     write ( *, '(a)' ) ' '
     write ( *, '(a,i8,a,a)' ) '  Process ', id, ' says "Hello, world!" ',name(1:clen)

     call MPI_Finalize ( error )
   end program

As an example, on a MAC with 2 cores that has mpich with gnu fortran you would issue the following two commands:
