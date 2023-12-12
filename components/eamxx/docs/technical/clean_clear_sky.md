# Clean- and clean-clear-sky diagnostics

In order to decompose the aerosol effective radiative forcing, additional diagnostic radiation calls are needed.
These extra diagnostics are optionally added to the main radiation call. The extra diagnostics are:

- Clean-clear-sky fluxes: the fluxes that would be present if there were neither aerosols nor clouds, and are calculated by adding an additional radiation call at the very beginning of the logic before the optics class is endowed with aerosol and cloud properties.
- Clean-sky fluxes: the fluxes that would be present if there were no aerosols, and are calculated by adding an additional radiation call after substantiating an additional optics class, but not endowing it with aerosol properties.

It was necessary to add an additional optics class because the original optics class is endowed with aerosols before clouds (in order to calculate the clear-sky fluxes).
The extra calls are controlled by runtime flags `extra_clnclrsky_diag` and `extra_clnsky_diag` (they take either `true` or `false` as their values).
