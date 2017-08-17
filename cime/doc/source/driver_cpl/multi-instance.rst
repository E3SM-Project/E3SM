Multi-instance Functionality
=============================

The multi-instance feature allows multiple instances of a given component to run in a single CESM run.
This might be useful for data assimilation or to average results from multiple instances to force another model.

The multi-instance implementation is fairly basic at this point.
It does not do any averaging or other statistics between multiple instances, and it requires that all prognostic components must run the same multiple instances to ensure correct coupling. 
The multi-instance feature is set via the ``$CASEROOT/env_mach_pes.xml`` variables that have an ``NINST_`` prefix.
The tasks and threads that are specified in multi-instance cases are distributed evenly between the multiple instances.
In other words, if 16 tasks are requested for each of two atmosphere instances, each instance will run on 8 of those tasks. 
The ``NINST_*_LAYOUT`` value should always be set to *concurrent* at this time.
Sequential running on multiple instances is still not a robust feature.
Multiple instances is a build time setting in env_mach_pes.xml. 
Multiple instance capabilities are expected to be extended in the future.
