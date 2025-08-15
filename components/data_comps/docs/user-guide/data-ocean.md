# The E3SM Data Ocean Model

<!-- disable certain linter checks here for more readable nested markdown  -->
<!-- markdownlint-disable  MD007 --> <!-- ul-indent -->

The E3SM data ocean has several different modes to support various realistic and idealized experiments. Sea surface temperatures (SST) can be either prescribed or prognostic. Prescribed SSTs are specified either through a data stream or analytically. Prognostic modes allow the SST field to evolve and respond to atmospheric conditions. The guides below provide more details on how to use these capabilities.

- Prescribed
    - [SST from Observations](#sst-from-observations)
    - [Idealized SST](#idealized-sst)
- Prognostic
    - [Traditional Slab Ocean Model](#traditional-slab-ocean-model) (SOM)
    - [Relaxed Slab Ocean](#relaxed-slab-ocean) (RSO)

## SST from Observations

Using SST data derived from observations is the most common use of the data ocean model, often for AMIP style experiments to reproduce historical periods.

Example compsets that use this capability are `F2010` and `F20TR`. These compsets use the `_DOCN%DOM_` compset modifier, which sets the `DOCN_MODE` variable in `env_run.xml`  to "prescribed".

Several additional XML variables need to be set in order to use this capability, which are set to defaults for common configurations, such as `F2010` at `ne30pg2` atmospheric resolution.

```text
SSTICE_DATA_FILENAME  Prescribed SST and ice coverage data file name
SSTICE_GRID_FILENAME  Grid file in "domain" format corresponding to SSTICE_DATA_FILENAME
SSTICE_YEAR_ALIGN     The model year that corresponds to SSTICE_YEAR_START on the data file
SSTICE_YEAR_START     The first year of data to use from SSTICE_DATA_FILENAME
SSTICE_YEAR_END       The last year of data to use from SSTICE_DATA_FILENAME
```

Most users will not need to edit these values from their defaults, but many scenarios require non-standard SST data, such as tropical cyclone hindcasts where the daily evolution of high-resolution SST data may be desireable.

## Idealized SST

The two main uses of idealized SST modes are aquaplanet (AQP) and radiative-convective equilibrium (RCE). The latter is just a special case of an aquaplanet where the SST is [usually] a constant value everywhere, traditionally used in conjunction with special modifications to homogenize radiation and disable rotation. There are several analytically specified SST patterns established by model intercomparison projects such as the Aqua-Planet Experiment (APE)[@blackburn_APE_context_2013] and RCEMIP[@wing_rcemip1_2018,@wing_rcemip2_2024].

### Idealized SST compsets

The following list shows the currently defined E3SM compsets that utilize idealized SST.

```text
FAQP
FAQP-MMF1
FAQP-MMF2
F-SCREAM-LR-AQP1
F-SCREAM-HR-AQP1
FRCE
FRCE-MMF1
FRCE-MMF2
FRCE-MW_295dT1p25
FRCE-MW_300dT0p625
FRCE-MW_300dT1p25
FRCE-MW_300dT2p5
FRCE-MW_305dT1p25
FRCE-MW-MMF1_295dT1p25
FRCE-MW-MMF1_300dT0p625
FRCE-MW-MMF1_300dT1p25
FRCE-MW-MMF1_300dT2p5
FRCE-MW-MMF1_305dT1p25
```

These all use "analytic" SST patterns that are specified via the `docn_comp_run()` subroutine in `components/data_comps/docn/src/docn_comp_mod.F90`. The `AQP` compsets currently only use the basic aquaplanet pattern that is symmetric about the equator. Other APE patterns introduce different meridional gradients and/or asymmetries. The various analytic SST patterns can be selected by changing the data ocean specifier: `_DOCN%AQP1_`.

The first 10 analytic aquaplanet SST patterns correspond to the aqua-planet experiment (APE) protocol as follows

```text
AQP1    = control symmetric SST pattern
AQP2    = Flat
AQP3    = Qobs = average of AQP1 and AQP2
AQP4    = Peaked
AQP5    = Control+5N
AQP6    = 1KEQ - small warm pool
AQP7    = 3KEQ - small warm pool
AQP8    = 3KW1 - large warm pool
AQP9    = Control+10N
AQP10   = Control+15N
```

!!!NOTE
    When using aquaplanet mode the orbital parameters will take on the idealized values shown below such that there are no seasonal variations, but there is still a diurnal cycle.
    ```text
    orb_eccen = 0
    orb_obliq = 0
    orb_mvelp = 0
    orb_mode = "fixed_parameters"
    ```

The basic RCE compsets use the `_DOCN%AQPCONST_` modifier to produce a globally constant SST value, which is set by the `DOCN_AQPCONST_VALUE` variable in `env_run.xml`. The "FRCE-MW" compsets were designed for RCEMIP-II to produce a "mock walker-cell" configuration, in which sinusoidal SST variations are applied to promote a coherent large-scale circulation.

### SST Data File

In addition to the analytic SST modes the user can also specify an idealized aquaplanet SST pattern via the `_DOCN%AQPFILE_` option. The `aquapfile` namelist variable is used to specify the SST pattern in this mode. Note that this option has not been used or tested recently, so the user may experience difficulty trying to use this feature.

## Traditional Slab Ocean Model

A slab ocean model (SOM) allows responsive SSTs to address the "infinite heat source" problem associated with prescribed SSTs, but is much cheaper than running with a full ocean model. The traditional SOM appraoch requires special inputs, such as a specified mixed layer depth pattern that can vary in time and a prescribed heat flux to account for the missing effects of ocean dynamics often referred to as "Q-flux". The Q-flux data is often estimated from a fully coupled simulation with active ocean and sea-ice so that the SOM simulation will resemble the full model.

Currently, we do not have Q-flux data to drive the SOM in E3SM. An alternative appraoch is to use a "relaxed" slab ocean (RSO) in which a specified relaxation time scale is used to bring the SST field back to a target SST field. The RSO mode is much simpler to use, but carries caveats that the user should be aware of before using. See [Data Ocean - Relaxed Slab Ocean](#relaxed-slab-ocean) for more information.

## Relaxed Slab Ocean

The relaxed slab ocean (RSO) is similar in many ways to the [traditional slab ocean model](#traditional-slab-ocean-model), but uses a specified relaxation time scale to avoid the need for specified "Q-flux" data to represent the effects of ocean transport. The RSO implementation in E3SM was inspired by Zarzycki (2016)[@Zarzycki_TC-ocn-cpl_2016].

A key consideration for the user is whether they need to use a realistic distribution of mixed layer depths (MLD), or whether their use case can benefit from the simplicity of a globally uniform MLD.

The RSO mode has the following namelist variables to influence the ocean behavior:

```text
RSO_relax_tau   SST relaxation timescale
RSO_fixed_MLD   globally uniform MLD value (use -1 for realistic MLD)
```

Other RSO parameter values are hardcoded in `components/data_comps/docn/src/docn_comp_mod.F90`.

```text
RSO_slab_option = 0           ! Option for setting RSO_X_cool
RSO_R_cool      = 11.75/86400 ! base cooling rate [K/s]
RSO_Tdeep       = 271.00      ! deep water temperature [K]
RSO_dT_o        = 27.0        ! scaling temperature gradient
RSO_h_o         = 30.0        ! scaling mixed layer depth
```

The RSO mode uses the `SSTICE_DATA_FILENAME` in `env_run.xml` for its data stream. For a globally uniform MLD this file only need to contain a `SST_cpl` variable for the SST that will act as the target SST value for relaxation. If a realistic MLD pattern is desired then the `hblt` variable must also be present. This data can be derived a number of ways, but we currently do not have a dedicated tool or workflow.
