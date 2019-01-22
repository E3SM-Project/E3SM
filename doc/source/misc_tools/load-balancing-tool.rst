.. _load_balancing_tool:


=========================
 CIME Load Balancing Tool
=========================

 Originally Developed by Sheri Mickelson mickelso@ucar.edu
 and Yuri Alekseev (ALCF/Argonne National Laboratory

 Updated 2017 Jason Sarich sarich@mcs.anl.gov (Argonne National Laboratory)


This Load Balancing tool performs several operations intended to find
a reasonable PE layout for CIME simulations. These operations involve two
steps::

  1. load_balancing_submit.py
     Run a series of simulations in order to obtain timing data

  2. load_balancing_solve.py
     Using the data provided in the previous program, solve a mixed integer
     linear program to optimize the model throughput. Requires installation
     of PuLP and uses the included COIN-CBC solver. (https://pythonhosted.org/PuLP)

Also in this documentation is::

  3. More about the algorithm used

  4. Extending the solver for other models

  5. Testing information for developers


*For the impatient*


1. set PYTHONPATH to include $CIME_DIR/scripts:$CIME_DIR/tools/load_balancing_tool
   
2. create PE XML <PESFILE> file to describe the PE layouts for the timing runs

3. $ ./load_balancing_submit.py --res <RESOLUTION> --compset <COMPSET> --pesfile <PESFILE>
   
4. ...  wait for jobs to run ...
   
5. $ ./load_balancing_solve.py --total-tasks <N> --blocksize 8



******************************************************************
Running simulations using load_balancing_submit.py
******************************************************************

Simulations can be run on a given system by executing the load_balancing_tool.py
script, located in cime/tools/load_balancing_tool/load_balancing_tool_submit.py.
This creates timing files in the case directory which will be used to solve
a mixed integer linear program optimizing the layout. If there is already timing
information available, then a 

As with the create_newcase and create_test scripts, command line options
are used to tailor the simulations for a given model. These values will be
directly forwarded to the passed::

     --compiler
     --project
     --compset (required)
     --res     (required)
     --machine

Other options include::

     --pesfile <filename>    (required)
          This file is used to designated the pes layout that
	  are used to create the timing data. The format is the same used
	  by CIME pes_files, but note that the 'pesize' tag will be used
	  to generate the casename. Also, this file will not be directly
	  passed through to CIME, but rather it will trigger xmlchange
	  commands to execute based on the values in the file.

     --test-id <prefix>
          By default, the load balancing tool will use casenames:
	     PFS_I0.res.compset.lbt
	     PFS_I1.res.compset.lbt
	     ...
	     PFS_IN.res.compset.lbt
	  for each simulation requested. These casenames will be forwarded to
	  the create_test script.

	  Using this option will instead direct the tool to use:
	     PFS_I0.res.compset.test-id
	     PFS_I1.res.compset.test-id
	     ...
	     PFS_IN.res.compset.test-id

     --force-purge
          Force the tool to remove any existing case directories if they
	  exist. Removes PFS_I*.res.compset.test-id

     --extra-options-file
          Add extra xml options to the timing runs from a user file,
	  these options will be set after create_newcase and before
	  case.setup.
	  This text file should have one variable per line in
	  the format <var>=<value>. Example:

	  STOP_OPTION=ndays
	  STOP_N=7
	  DOUT_S=FALSE


******************************************************************
Optimizing the layout using load_balacing_solve.py
******************************************************************

Reads timing data created with load_balancing_submit.py (or otherwise,
see --timing-files option) and solves an mixed integer optimization problem
using these timings. The default layout (IceLndAtmOcn) minimizes the cost per
model day assuming the layout::

              ____________________
             | ICE  |  LND  |     |
             |______|_______|     |
             |              | OCN |
             |    ATM       |     |
             |______________|_____|


An IceLndWavAtmOcn layout is also available.  It is possible to extend
this tool to solve for other layouts (See Section 1.4 Extending the Load
Balancing Tool)

Note -- threading is not considered part of this optimization, it is assumed that
all timing data have the same threading structure (i.e. all ATM runs use two threads per PE)

Options recognized by the solver::

  --layout <class_name>
      Name of the class used to solve the layout problem. The only built-in
      class at this time is the default IceLndAtmOcn, but this can be extended.
      See section 4 Extending the Load Balancing Tool

  --total-tasks N    (required)
      The total number of PEs that can be assigned

  --timing-dir <dir>
      Optional, read in all files from this directory as timing data

  --test-id <prefix>
      The test-id used when submitting the timing jobs. This option can also
      be used to set a single directory where ALL of the timing data is.
      The solver will extract data from timing files that match either pattern:
         <prefix>.test-id/timing/timing.<prefix>.test-id
	 <prefix>.test-id/timing/timing.<prefix>.test-id

  --blocksize N
      The blocksize is the granularity of processors that will be group
      together, useful for when PEs to be multiples of 8, 16, etc.

  --blocksize-XXX N
      Components don't all have to have the same blocksize. The default
      blocksize given by --blocksize can be overridden for a given component
      using this option, where XXX can be ATM, ICE, GLC, etc.
      Example:
      --blocksize 8 --blocksize-GLC 1
          will set the GLC blocksize to 1 and all other blocksizes to 8

  --milp-output <filename>
      After extracting data from timing files and before solving, write the
      data to a .json file where is can be analyzed or manually edited.

  --milp-input <filename>
      Read in the problem from the given .json file instead of extracting from
      timing files.

  --pe-output <filename>
      Write the solution PE layout to a potential pe xml file.


***************************
More about the algorithm
***************************

Before solving the mixed-integer linear program, a model of the cost vs ntasks
function is constructed for each component.

Given a component data set of costs (C1,C2,..,Cn) and nblocks (N1,N2,..,Nn),
then an piecewise set of n+1 linear constraints are created using the idea:

If N < N1 (which means that N1 cannot be 1), then assume that there is
perfect scalability from N to N1. Thus the cost is on the line
defined by the points (1, C1*N1) - (N1, C1).

If N is between N_i and N_{i+1}, then the cost is on the line defined by the
points (N_i, C_i) and (N_{i+1}, C_{i+1}.

If N > Nn, then we want to extrapolate the cost at N=total_tasks
(we define N{n+1} = total_tasks, C{n+1} = estimated cost using all nodes)
Assuming perfect scalability is problematic at this level, so we instead
assume that the parallel efficiency drops at the same factor as it does

  from N=N{n-1} to N = Nn

      First solve for efficiency E:
      C{n-1} - Cn = E * (C{n-1} * N{n-1} / Nn)

      Then E to find C{n+1} (cost at ntasks N{n+1}):
      Cn - Ct = E * (Cn * Nn / Nt)

      Now cost is on the line defined by (Nn,Cn) - (Nt,Ct)

Assuming that this piecewise linear function describes a convex function, we do
not have to explicitly construct this piecewise function and can instead use
each of the cost functions on the entire domain.

These piecewise linear models give us the following linear constraints, where
the model time cost C as a function of N (ntasks) for each component
is constrained by::

  C >= Ci  - Ni * (C{i+1}-Ci) / (N{i+1}-Ni) +
             N *  (C{i+1}-Ci) / (N{i+1}-Ni)    for i=0..n


These constraints should be in effect for any extensions of the solver (the
components involved may be different).

There are options available in load_balancing_submit.py to inspect these
piecewise linear models::

	  --graph-models (requires matplotlib)
	  --print-models (debugging modes writes the models to the log)


Now that these constraints are defined, the mixed integer linear program (MILP)
follows from the layout::

     NOTES: variable N[c] is number of tasks assigned for component c
            variable NB[c] is the number of blocks assigned to component c
            constant C[c]_i is the cost contributed by component c from
	                  timing data set i
            constant N[c]_i is the ntasks assigned to component c from
	                  timing data set i

              ____________________
             | ICE  |  LND  |     |
       T1    |______|_______|     |
             |              | OCN |
             |    ATM       |     |
       T     |______________|_____|

      Min T
      s.t.  Tice      <= T1
            Tlnd      <= T1
            T1 + Tatm <= T
            Tocn      <= T

            NB[c]        >= 1 for c in [ice,lnd,ocn,atm]
            N[ice] + N[lnd] <= N[atm]
            N[atm] + N[ocn] <= TotalTasks
	    N[c] = blocksize * NB[c], for c in [ice,lnd,ocn,atm]


            T[c]        >= C[c]_{i} - N[c]_{i} *
                       (C[c]_{i+1} - C[c]_{i}) / (N[c]_{i+1} - N[c]_{i})
                       + N[c] * (C[c]_{i+1} - C[c]_{i})
                                               / (N[c]_{i+1} - N[c]_{i}),
                        for i=0..#data points (original + extrapolated,
		            c in [ice,lnd,ocn,atm]
            all T vars >=0
	    all N,NB vars integer

This MILP is solved using the PuLP python interface to the COIN-CBC solver
https://pythonhosted.org/PuLP/
https://www.coin-or.org/Cbc/


************************************
Extending the Load Balancing Tool
************************************
The file $CIME_DIR/tools/load_balancing_tool/optimize_model.py
contains a base class OptimizeModel as well as an implementation class
IceLndAtmOcn. Any layout solver will look similar to IceLndAtmOcn
except for the components involved and the layout-specific constraints.

Example class and inherited methods that should be overridden:

file my_new_layout.py::

  import optimize_model

  class MyNewLayout(optimize_model.OptimizeModel)
     def get_required_components(self):
         """
         Should be overridden by derived class. Return a list of required
         components (capitalized) used in the layout.
         Example: return ['ATM', 'LND', 'ICE']
         """

     def optimize(self):
          """
          Run the optimization.
          Must set self.state using LpStatus object
          LpStatusOptimal    -> STATE_SOLVED_OK
          LpStatusNotSolved  -> STATE_UNSOLVED
          LpStatusInfeasible -> STATE_SOLVED_BAD
          LpStatusUnbounded  -> STATE_SOLVED_BAD
          LpStatusUndefined  -> STATE_UNDEFINED
          -- use self.set_state(lpstatus) --
          Returns state

          If solved, then solution will be stored in self.X dictionary, indexed
          by variable name. Suggested convention:
          'Tice', 'Tlnd', ... for cost per component
          'Nice', 'Nlnd', ... for ntasks per component
          'NBice', 'NBlnd', ... for number of blocks per component

          The default implementation of get_solution() returns a dictionary
          of these variable keys and their values.
          """

     def get_solution(self):
         """
         Return a dictionary of the solution variables, can be overridden.
         Default implementation returns values in self.X
         """


To use this new layout:
   1. save the class MyNewLayout in file my_new_layout.py
   2. make sure that my_new_layout.py is in PYTHONPATH
   3. Use those names in your execution command line argument to --layout
      ::

         $ ./load_balancing_solve.py ... --layout my_new_layout.MyNewLayout

To permanently add to CIME:

   1. add MyNewLayout class to layouts.py
   2. run using '--layout MyNewLayout'
   3. add test in tests/load_balance_test.py that uses that name in command
      line argument (see test for atm_lnd)
   4. make pull request


*******
Testing
*******

To run the provided test suite: 

  1. set PYTHONPATH to include CIME libraries::

      $ export CIME_DIR=/path/to/cime
      $ export PYTHONPATH=$CIME_DIR/scripts:$CIME_DIR/tools/load_balancing_tool

  2. To run an example::

      $ cd $CIME_DIR/tools/load_balancing_tool
      $ ./load_balancing_solve.py --json-input tests/example.json --blocksize 8
      Solving Mixed Integer Linear Program using PuLP interface to COIN-CBC
      PuLP solver status: Solved
      COST_ATM = 22.567587
      COST_ICE = 1.375768
      COST_LND = 1.316000
      COST_OCN = 15.745000
      COST_TOTAL = 23.943355
      NBLOCKS_ATM = 124
      NBLOCKS_ICE = 109
      NBLOCKS_LND = 15
      NBLOCKS_OCN = 4
      NTASKS_ATM = 992
      NTASKS_ICE = 872
      NTASKS_LND = 120
      NTASKS_OCN = 32
      NTASKS_TOTAL = 1024

  3. To run the test suite::

	$ cd $CIME_DIR/tools/load_balancing_tool
        $ ./tests/load_balancing_test.py


      
