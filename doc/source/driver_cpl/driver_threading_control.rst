Driver Threading Control
========================

OpenMP thread counts are controlled at three levels.

- The coarsest level is prior to launching the model. The environment variable OMP_NUM_THREADS is usually set to the largest value any mpi task will use. At a minimum, this will ensure threading is turned on to the maximum desired value in the run. 

- The next level is during the driver initialization phase. When the mpi communicators are initialized, the maximum number of threads per mpi task can be computed based on the ccsm_pes namelist input.  At that point, there is an initial fortran call to the intrinsic, omp_set_num_threads. When that happens and if that call is successful, the number of threads will be set to the maximum needed in the system on an mpi task by task basis. 

- Finally, there is the ability of CESM to change the thread count per task as each component is individually called and as the model integrates through the driver run loop. In other words, for components that share the same hardware processor but have different threads per task, this feature allows those components to run with the exact value set by the user in the ccsm_pes namelist.  This final level of thread control is turned off by default, but it can be turned on using the driver namelist variable ``drv_threading``.

This fine control of threading is likely of limited use at this point given the current driver implementation.
