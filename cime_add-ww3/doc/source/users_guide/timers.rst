.. _timers:

===================
Timers and timing
===================

CIME includes a copy of the General Purpose Timing Library (GPTL) and timers are placed throughout the CIME driver.  CIME-driven models typically
also have GPTL timers in their code and very detailed timing information can be obtained.

.. _model-timing-data:

Model timing data
------------------

Every model run produces three types of timing output that you can examine:

1. **$CASEROOT/timing/$model_timing.$CASE.$datestamp**

   This is the most useful way to quickly determine timing summaries across components
   The following describes the most important parts of this timing file:

   An example timing file of this type is:

   ::

      ---------------- TIMING PROFILE ---------------------
      Case        : b.e20.BHIST.f09_g17.20thC.297_02
      LID         : 9459679.chadmin1.180517-114852
      Machine     : cheyenne
      Caseroot    : /glade/p/cesmdata/cseg/runs/cesm2_0/b.e20.BHIST.f09_g17.20thC.297_02
      Timeroot    : /glade/p/cesmdata/cseg/runs/cesm2_0/b.e20.BHIST.f09_g17.20thC.297_02/Tools
      User        : hannay
      Curr Date   : Thu May 17 12:42:27 2018
      grid        : a%0.9x1.25_l%0.9x1.25_oi%gx1v7_r%r05_g%gland4_w%ww3a_m%gx1v7
      compset     : HIST_CAM60_CLM50%BGC-CROP_CICE_POP2%ECO_MOSART_CISM2%NOEVOLVE_WW3_BGC%BDRD
      run_type    : hybrid, continue_run = FALSE (inittype = TRUE)
      stop_option : nyears, stop_n = 1
      run_length  : 365 days (364.958333333 for ocean)

      component       comp_pes    root_pe   tasks  x threads instances (stride)
      ---------        ------     -------   ------   ------  ---------  ------
      cpl = cpl        3456        0        1152   x 3       1      (1     )
      atm = cam        3456        0        1152   x 3       1      (1     )
      lnd = clm        2592        0        864    x 3       1      (1     )
      ice = cice       864         864      288    x 3       1      (1     )
      ocn = pop        768         1152     256    x 3       1      (1     )
      rof = mosart     2592        0        864    x 3       1      (1     )
      glc = cism       3456        0        1152   x 3       1      (1     )
      wav = ww         96          1408     32     x 3       1      (1     )
      esp = sesp       1           0        1      x 1       1      (1     )

      total pes active           : 12960
      mpi tasks per node               : 36
      pe count for cost estimate : 4320

      Overall Metrics:
      Model Cost:            3541.30   pe-hrs/simulated_year
      Model Throughput:        29.28   simulated_years/day

      Init Time   :     242.045 seconds
      Run Time    :    2951.082 seconds        8.085 seconds/day
      Final Time  :       0.008 seconds

      Actual Ocn Init Wait Time     :     768.737 seconds
      Estimated Ocn Init Run Time   :       0.248 seconds
      Estimated Run Time Correction :       0.000 seconds
        (This correction has been applied to the ocean and total run times)

      Runs Time in total seconds, seconds/model-day, and model-years/wall-day
      CPL Run Time represents time in CPL pes alone, not including time associated with data exchange with other components

      TOT Run Time:    2951.082 seconds        8.085 seconds/mday        29.28 myears/wday
      CPL Run Time:     248.696 seconds        0.681 seconds/mday       347.41 myears/wday
      ATM Run Time:    2097.788 seconds        5.747 seconds/mday        41.19 myears/wday
      LND Run Time:     545.991 seconds        1.496 seconds/mday       158.24 myears/wday
      ICE Run Time:     389.173 seconds        1.066 seconds/mday       222.01 myears/wday
      OCN Run Time:    2169.399 seconds        5.944 seconds/mday        39.83 myears/wday
      ROF Run Time:      42.241 seconds        0.116 seconds/mday      2045.41 myears/wday
      GLC Run Time:       1.049 seconds        0.003 seconds/mday     82364.16 myears/wday
      WAV Run Time:     517.414 seconds        1.418 seconds/mday       166.98 myears/wday
      ESP Run Time:       0.000 seconds        0.000 seconds/mday         0.00 myears/wday
      CPL COMM Time:   2464.660 seconds        6.752 seconds/mday        35.06 myears/wday

      ---------------- DRIVER TIMING FLOWCHART ---------------------
      .............


   TIMING PROFILE is the first section in the timing output. It
   summarizes general timing information for the run. The total run
   time and cost are given in several metrics to facilitate analysis
   and comparisons with other runs. These metrics includ pe-hrs per
   simulated year (cost), simulated years per wall day (thoughput),
   seconds, and seconds per model day. The total run time for each
   component and the time for initialization of the model also are
   provided. These times are the aggregate over the total run and do
   not take into account any temporal or processor load imbalances.

   DRIVER TIMING FLOWCHART is the second section in the timing
   output. It provides timing information for the driver in
   sequential order and indicates which processors are involved in
   the cost. Finally, the timings for the coupler are broken out at
   the bottom of the timing output file.


2. **$CASEROOT/timing/$model_timing_stats.$date**

   Provides an overall detailed timing summary for each component, including the minimum and maximum of all the model timers.

3. **cpl.log.$datestamp**

   Contains the run time for each model day during the run and is
   output during the run. You can search for ``tStamp`` in the cpl.log
   file to see the information, which is useful for tracking down
   temporal variability in cost due to inherent model variability or
   to hardware. The model daily cost generally is pretty constant
   unless I/O is written intermittently, such as at the end of the
   month. This file will appear either in **$RUNDIR** or in
   **DOUT_S_ROOT/logs** for your run.

The xml variable ``CHECK_TIMING``, if set to ``TRUE`` (the default) will produce the timing files in the **$CASEROOT/timing** directory.


Controlling timers
------------------

User customization of timers is done via the xml variables ``TIMER_LEVEL`` and ``TIMER_DETAIL``.

* ``TIMER_LEVEL``:

  This is the maximum code stack depth of enabled timers.

* ``TIMER_DETAIL``:

  This is an integer indicating maximum detail level to profile.  This
  xml variable is used to set the namelist variable timing_detail_limit.
  This namelist variable is used by perf_mod (in
  $CIMEROOT/src/share/timing/perf_mod.F90) to turn timers off and on
  depending on calls to the routine t_adj_detailf.  If in the code a
  statement appears like t_adj_detailf(+1), then the current timer
  detail level is incremented by 1 and compared to the time_detail_limit
  obtained from the namelist.  If the limit is exceeded then the timer
  is turned off.

Further control of timers is then done via modifications of the **prof_inparm namelists** in the file **drv_in**. This is done
via keyword-value settings in user_nl_cpl. As an example, if you want to set the namelist variable ``profile_barriers`` to ``.true.``,
add the following line in your **$CASEROOT/user_nl_cpl**:

::

   profile_barriers = .true.


Advice on setting your wallclock time
-------------------------------------

When you look at the **$model_timing.$CASE.$datestamp** file for "Model Throughput", you will find output like this:
 ::

  Overall Metrics:
  Model Cost: 327.14 pe-hrs/simulated_year (scale= 0.50)
  Model Throughput: 4.70 simulated_years/day

The model throughput is the estimated number of model years that you
can run in a wallclock day. Based on this, you can maximize your queue
limit and change ``$STOP_OPTION`` and ``$STOP_N``.

For example, say a model's throughput is 4.7 simulated_years/day, and
the maximum runtime limit on your machine is 12 hours. 4.7 model
years/24 hours * 12 hours = 2.35 years. On the massively parallel
computers, there is always some variability in how long it will take
a job to run. On some machines, you may need to leave as much as 20%
buffer time in your run to guarantee that jobs finish reliably before
the time limit. For that reason, set your model to run only one model
year/job. In this example, set your wallclock at 12 hours and invoke
`xmlchange <../Tools_user/xmlchange.html>`_  in ``CASEROOT`` as shown here: ::

  >./xmlchange STOP_OPTION=nyears
  >./xmlchange STOP_N=1
  >./xmlchange REST_OPTION=nyears
  >./xmlchange REST_N=1
