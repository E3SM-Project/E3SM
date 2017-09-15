History and Restarts
====================

In addition to log files, component models also produce history and restart files.
History files are generally netcdf format and contain fields associated with the state of the model. 
History files are implemented and controlled independently in the component models, although support for monthly average history files is a standard output of most production runs.
The driver has a file naming standard for history files which includes the case names, component name, and model date.

All CIME-compliant component models must be able to stop in the middle of a run and then subsequently restart in a bit-for-bit fashion.
For most models, this requires the writing of a restart file. 
The restart file can be any format, although netcdf has become relatively standard, and it should contain any scalars, fields, or information that is required to restart the component model in exactly the same state as when the restart was written and the model was stopped.
The expectation in CIME is that a restart of a model run will be bit-for-bit identical and this is regularly tested as part of component model development by running the model 10 days, writing a restart at the end of 5 days, and then restarting at day 5 and comparing the result with the 10 day run. 
Unlike history files, restart files must be coordinated across different components.
The restart frequency is set in the driver time manager namelist by driver namelist variables ``restart_option``, ``restart_n``, and ``restart_ymd``.
The driver will trigger a restart alarm in clocks when a coordinated restart is requested. 
The components are required to check this alarm whenever they are called and to write a restart file at the end of the current coupling period.
This method ensures all components are writing restart files at a consistent timestamp. 
The restart filenames are normally set in a generic rpointer file.
The rpointer file evolves over the integration and keeps track of the current restart filenames. 
When a model is restarted, both the rpointer file and the actual restart file are generally required.

Many models are also able to restart accumulating history files in the middle of an accumulation period, but this is not a current requirement for CIME compliant components.
In production, the model is usually started and stopped on monthly boundaries so monthly average history files are produced cleanly. 
The run length of a CESM1 production run is usually specified using the nmonths or nyears option and restart files are normally written only at the end of the run.
