# Data Ocean - SST from Observations

Using SST data derived from observations is the most common use of the data ocean model, often for AMIP style experiments to reproduce historical periods.

Example compsets that use this capability are `F2010` and `F20TR`. These compsets use the `_DOCN%DOM_` compset modifier, which sets the `DOCN_MODE` variable in `env_run.xml`  to "prescribed".

Several additional XML variables need to be set in order to use this capability, which are set to defaults for common configurations, such as `F2010` at `ne30pg2` atmospheric resolution.

```
SSTICE_DATA_FILENAME  Prescribed SST and ice coverage data file name
SSTICE_GRID_FILENAME  Grid file in "domain" format corresponding to SSTICE_DATA_FILENAME
SSTICE_YEAR_ALIGN     The model year that corresponds to SSTICE_YEAR_START on the data file
SSTICE_YEAR_START     The first year of data to use from SSTICE_DATA_FILENAME
SSTICE_YEAR_END       The last year of data to use from SSTICE_DATA_FILENAME
```

Most users will not need to edit these values from their defaults, but many scenarios require non-standard SST data, such as tropical cyclone hindcasts where the daily evolution of high-resolution SST data may be desireable.