# Atmospheric Processes

In EAMxx, the `AtmosphereProcess` (AP) is a class representing a portion of the atmosphere timestep algorithm.
In simple terms, an AP is an object that given certain input fields performs some calculations to compute
some output fields.

TODO: describe init sequcene (e.g., the process of requesting fields), base class main
      interfaces/capabilities (e.g., subcycling), class expectations (e.g., must update fields on physics grid)

Here is a list of currently implemented atmosphere processes.
TODO: add links to papers/github-repos, and a SMALL description
* p3: Microphysics, blah blah
* SHOC: Macrophysics/Turbulence, blah
* rrtmgp: Radiation, blah
* spa: prescribed aerosols, blah blah
* surface coupling: blah
* mam: prognostic aerosols, blah blah
* nudging: This process is responsible for nudging the model simulation given a set of files with a target nudged state.
