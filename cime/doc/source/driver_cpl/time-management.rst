.. _time-management:

===============
Time Management
===============

-------------
Driver Clocks
-------------
The driver manages the main clock in the system.  The main clock
advances at the shortest coupling period and uses alarms to trigger
component coupling and other events.  In addition, the driver
maintains a clock that is associated with each component.  The
driver's component clocks have a timestep associated with the coupling
period of that component.  The main driver clock and the component
clocks in the driver advance in a coordinated manor and are always
synchronized.  The advancement of time is managed as follows in the
main run loop.  First, the main driver clock advances one timestep and
the component clocks are advanced in a synchronous fashion.  The clock
time represents the time at the end of the next model timestep.
Alarms may be triggered at that timestep to call the the atmosphere,
land, runoff, sea ice, land ice, or ocean run methods.  If a component
run alarm is triggered, the run method is called and the driver passes
that component's clock to that component.  The component clock
contains information about the length of the next component
integration and the expected time of the component at the end of the
integration period.

Generally, the component models have indepedent time management
software.  When a component run method is called, the component must
advance the proper period and also check that their internal clock is
consistent with the coupling clock before returning to the driver.
The clock passed to the component by the driver contains this
information.  Component models are also responsible for making sure
the coupling period is consistent with their internal timestep.
History files are managed independently by each component, but restart
files are coordinated by the driver.

The driver clocks are based on ESMF clock datatype are are supported
in software by either an official ESMF library or by software included
in CIME called ``esmf_wrf_timemgr``, which is a much simplified
Fortran implementation of a subset of the ESMF time manager
interfaces.

--------------------
The Driver Time Loop
--------------------
The driver time loop is hardwired to sequence the component models in
a specific way to meet scientific requirements and to otherwise
provide the maximum amount of potential concurrency of work.  The
results of the model integration are not dependent on the processor
layout of the components.  See the Craig et al IJHPCA 2012 reference
for further details.

In addition, the driver is currently configured to couple the
atmosphere, land, and sea ice models using the same coupling frequency
while the runoff, land ice, and ocean model can be coupled at the same
or at a lower frequency.  To support this feature, the driver does
temporal averaging of coupling inputs to the ocean and runoff, and the
driver also computes the surface ocean albedo at the higher coupling
frequency.  There is no averaging of coupling fields for other
component coupling interactions and the land and sea ice models'
surface albedos are computed inside those components.  Averaging
functionality could be added to the driver to support alternative
relative coupling schemes in the future if desired with the additional
caveat that the interaction between the surface albedo computation in
each component and the atmospheric radiation calculation have to be
carefully considered.  In addition, some other features may need to be
extended to support other coupling schemes and still allow model
concurrency.

The coupler processors (pes) handle the interaction of data between
components, so there are separate tasks associated with deriving
fields on the coupler pes, transfering data to and from the coupler
pes and other components, and then running the component models on
their processors.  The driver time loop is basically sequenced as
follows,
::

   - The driver clock is advanced first and alarms set.
   - Input data for ocean, land, sea ice, and runoff is computed.
   - Ocean data is rearranged from the coupler to the ocean pes.
   - Land data is rearranged from the coupler to the land pes.
   - Ice data is rearranged from the coupler to the ice pes.
   - Runoff data is rearranged from the coupler to the ice pes.
   - The ice model is run.
   - The land model is run.
   - The runoff model is run.
   - The ocean model is run.
   - The ocean inputs are accumulated, and the atmosphere/ocean fluxes are
     computed on the coupler pes based on the results from the previous
     atmosphere and ocean coupled timestep.
   - Land data is rearranged from the land pes to the coupler pes.
   - Land ice input is computed.
   - Land ice data is rearranged from the coupler to the land ice pes.
   - River output (runoff) data is rearranged from the runoff pes to the coupler pes.
   - Ice data is rearranged from the ice pes to the coupler pes.
   - Coupler fractions are updated.
   - Atmospheric forcing data is computed on the coupler pes.
   - Atmospheric data is rearranged from the coupler pes to the atmosphere pes.
   - The atmosphere model is run.
   - The land ice model is run.
   - Land ice data is rearranged from the land ice pes to the coupler pes.
   - Atmospheric data is rearranged from the atmosphere pes to the coupler pes.
   - Ocean data is rearranged from the ocean pes to the coupler pes.
   - The loop returns

    Within this loop, as much as possible, coupler work associated 
    with mapping data, merging fields, diagnosing, applying area corrections, 
    and computing fluxes is overlapped with component work.

The land ice model interaction is slightly different. 
::

   - The land ice model is run on the land grid
   - Land model output is passed to the land ice model every land coupling period.
   - The driver accumluates this data, interpolates the data to the land ice grid, 
     and the land ice model advances the land ice model about once a year.

The runoff coupling should be coupled at a frequency between the land
coupling and ocean coupling frequencies. The runoff model runs at the
same time as the land and sea ice models when it runs.

The current driver sequencing has been developed over nearly two
decades, and it plays a critical role in conserving mass and heat,
minimizing lags, and providing stability in the system.  The above
description is consistent with the `concurrency limitations 
<http://www.cesm.ucar.edu/models/cesm2.0/cpl7/doc.new/x32.html#design_seq>`_.
Just to reiterate, the land, runoff, and sea ice models will always
run before the atmospheric model, and the coupler and ocean models are
able to run concurrently with all other components.  The coupling
between the atmosphere, land, sea ice, and atmosphere/ocean flux
computation incurs no lags but the coupling to the ocean state is
lagged by one ocean coupling period in the system.  `Mass and heat
<http://www.cesm.ucar.edu/models/cesm2.0/cpl7/doc.new/x168.html>`_
are conserved in the system with more description.


It is possible to reduce the ocean lag in the system.  A driver
namelist variable, ``ocean_tight_coupling``, moves the step where
ocean data is rearranged from the ocean pes to the coupler pes from
the end of the loop to before the atmosphere/ocean flux computation.
If ocean_tight_coupling is set to true, then the ocean lag is reduced
by one atmosphere coupling period, but the ability of the ocean model
to run concurrently with the atmosphere model is also reduced or
eliminated.  This flag is most useful when the ocean coupling
frequency matches the other components.

------------------
Coupling Frequency
------------------
In the current implementation, the coupling period must be identical
for the atmosphere, sea ice, and land components.  The ocean coupling
period can be the same or greater.  The runoff coupling period should
be between or the same as the land and ocean coupling period.  All
coupling periods must be multiple integers of the smallest coupling
period and will evenly divide the NCPL_BASE_PERIOD, typically one day,
set in env_run.xml.  The coupling periods are set using the NCPL env
variables in env_run.xml.

The coupling periods are set in the driver namelist for each component
via variables called something like atm_cpl_dt and atm_cpl_offset.
The units of these inputs are seconds.  The coupler template file
derives these values from CIME script variable names like ATM_NCPL
which is the coupling frequency per day.  The \*_cpl_dt input
specifies the coupling period in seconds and the \*_cpl_offset input
specifies the temporal offset of the coupling time relative to initial
time.  An example of an offset might be a component that couples every
six hours.  That would normally be on the 6th, 12th, 18th, and 24th
hour of every day.  An offset of 3600 seconds would change the
coupling to the 1st, 7th, 13th, and 19th hour of every day.  The
offsets cannot be larger than the coupling period and the sign of the
offsets is such that a positive offset shifts the alarm time forward
by that number of seconds.  The offsets are of limited use right now
because of the limitations of the relative coupling frequencies.

Offsets play an important role in supporting concurrency.  There is an
offset of the smallest coupling period automatically introduced in
every coupling run alarm for each component clock.  This is only
mentioned because it is an important but subtle point of the
implementation and changing the coupling offset could have an impact
on concurrency performance.  Without this explicit automatic offset,
the component run alarms would trigger at the end of the coupling
period.  This is fine for components that are running at the shortest
coupling period, but will limit the ability of models to run
concurrently for models that couple at longer periods.  What is really
required for concurrency is that the run alarm be triggered as early
as possible and that the data not be copied from that component to the
coupler pes until the coupling period has ended.  The detailed
implementation of this feature is documented in the seq_timemgr_mod.
90 file and the impact of it for the ocean coupling is implemented in
the ccsm_driver.F90 code via use of the ocnrun_alarm and ocnnext_alarm
variables.

