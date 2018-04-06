.. _porting:

==============================================
Porting and validating CIME on a new platform
==============================================

One of the first steps for many users is getting CIME-based models running on their local machine.
This section describes that process.

Required libraries/packages
---------------------------

The machine needs to have:

- a functioning MPI environment (unless you plan to run on a single core with the CIME mpi-serial library).
- build tools gmake and cmake,
- a netcdf library version 4.3 or newer built with the same compiler you will use for CIME.

A pnetcdf library is optional.

If you are using MPI, make sure you can run a basic MPI parallel program on your machine before you attempt a CIME port. You can use this :ref:`MPI example <mpi-example>` to check.

.. _mpi-example:

An MPI example
---------------

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




Steps for porting
---------------------------

Porting CIME involves several steps in which you create, at a minimum, a **config_machines.xml** file in your **$HOME/.cime** directory.
In addition, if you have a batch system, you will also need to add a **config_batch.xml** file to your **$HOME/.cime** directory.

All files in **$HOME/.cime/** are appended to the xml objects that are read into memory from the **$CIME/config/[model]**, where **[model]** is either ``e3sm`` or ``cesm``.

Follow these steps:

#. Create a **$HOME/.cime** directory

#. Create a **config** file in that directory. It should contain the following two lines:
   ::

      [main]
      CIME_MODEL=cesm

   or

   ::

      [main]
      CIME_MODEL=e3sm

#. Create a **config_machines.xml** file in the same directory.

   This file contains all the information you need to set in order to configure a new machine to be CIME-compliant.

   Use one of the templates described here to create the file.

   - If you are on a MAC and have the libraries installed with either MacPorts or HomeBrew, copy
     **$CIME/config/cesm/machines/userdefined_laptop_template/config_machines.xml** to
     **$HOME/.cime/config_machines.xml**.

   - Otherwise, copy **$CIME/config/xml_schemas/config_machines_template.xml** to
     **$HOME/.cime/config_machines.xml**.

   Fill in the contents of **$HOME/.cime/config_machines.xml** that are specific to your machine. For more details see :ref:`the config_machines.xml file <machinefile>`.

   Check to ensure that your **config_machines.xml** file conforms to the CIME schema definition by doing the following:
   ::

      xmllint --noout --schema $CIME/config/xml_schemas/config_machines.xsd $HOME/.cime/config_machines.xml

#. If you have compiler settings that are specific to your machine, create a **$HOME/.cime/config_compilers.xml** file.

   The default compiler settings are set in **$CIME/config/[model]/machines/config_compilers.xml**, where **[model]** can be either ``e3sm`` or ``cesm``.

   There is no template for **config_compilers.xml**.

#.  If you have a batch system, create a **$HOME/.cime/config_batch.xml** file.

   Out-of-the-box batch settings are set in **$CIME/config/[model]/machines/config_batch.xml**, where **[model]** can be either ``e3sm`` or ``cesm``.

#. Once you have defined a basic configuration for your machine in your **$HOME/.cime** xml files, run **scripts_regression_test.py** interactively from the **$CIME/scripts/tests** directory.
   This performs a number of basic unit tests starting from the simplest and working toward more complicated ones.

After running those steps correctly, you are ready to try a case at your target compset and resolution.
   Once you have successfully created the required xml files in your .cime directory and are satisfied with the results you can merge them into the default files in the **config/$CIME_MODEL/machines** directory.
   If you would like to make this machine definition available generally you may then issue a pull request to add your changes to the git repository.

Validating your port
---------------------------

The following port validation is recommended for any new machine.
Carrying out these steps does not guarantee the model is running properly in all cases nor that the model is scientifically valid on the new machine.

In addition to these tests, detailed validation should be carried out for any new production run.
That means verifying that model restarts are bit-for-bit identical with a baseline run, that the model is bit-for-bit reproducible when identical cases are run for several months, and that production cases are monitored carefully as they integrate forward to identify any potential problems as early as possible. Users are responsible for their own validation process, especially with respect to science validation.

These are the recommended steps for validating a port:

1. Verify functionality by performing these `functionality tests <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_:

::

   ERS_D.f19_g16.X
   ERS_D.T31_g37.A
   ERS_D.f19_g16.B1850CN
   ERI.ne30_g16.X
   ERI.T31_g37.A
   ERI.f19_g16.B1850CN
   ERS.ne30_ne30.F
   ERS.f19_g16.I
   ERS.T62_g16.C
   ERS.T62_g16.DTEST
   ERT.ne30_g16.B1850CN


2. Verify performance and scaling analysis.

   a. Create one or two `load-balanced <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ configurations to check into ``Machines/config_pes.xml`` for the new machine.

   b. Verify that performance and scaling are reasonable.

   c. Review timing summaries in ``$CASEROOT`` for load balance and throughput.

   d. Review coupler "daily" timing output for timing inconsistencies.
      As mentioned in `load balancing a case <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_, useful timing information is contained in a **cpl.log.$date** file that is produced for every run.
      The file contains the run time for each model day during the model run.
      This diagnostic is output as the model runs.
      Searc for ``tStamp`` in this file to see this information.
      The timing information is useful for tracking down temporal variability in model cost due to either inherent model variability cost (I/O, spin-up, seasonal, and so on) or hardware.
      The model daily cost generally is pretty constant unless I/O is written intermittently, such as at the end of the month.

3. Perform validation (both functional and scientific):

   a. Perform a new CIME validation test (**TODO: fill this in**)

   b. Follow the `CCSM4.0 CICE port-validation procedure <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

   c. Follow the `CCSM4.0 POP2 port-validation procedure <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

4. Perform two, one-year runs (using the expected load-balanced configuration) as separate job submissions and verify that atmosphere history files are BFB for the last month.
   Do this after some performance testing is complete; you can also combine this with the production test by running the first year as a single run and the second year as a multi-submission production run.
   This will test reproducibility, exact restart over the one-year timescale, and production capability all in one test.

5. Carry out a 20- to 30-year 1.9x2.5_gx1v6 resolution, B_1850_CN compset simulation and compare the results with the diagnostics plots for the 1.9x2.5_gx1v6 Pre-Industrial Control (see the `CCSM4.0 diagnostics <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_).
   Model output data for these runs will be available on the `Earth System Grid (ESG) <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ as well.




