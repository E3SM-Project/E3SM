# Data Ocean - Idealized

The two main uses of idealized SST modes are aquaplanet (AQP) and radiative-convective equilibrium (RCE). The latter is just a special case of an aquaplanet where the SST is [usually] a constant value everywhere, traditionally used in conjunction with special modifications to homogenize radiation and disable rotation. There are several analytically specified SST patterns established by model intercomparison projects such as the Aqua-Planet Experiment (APE)[@blackburn_APE_context_2013] and RCEMIP[@wing_rcemip1_2018,@wing_rcemip2_2024].

## Idealized SST compsets

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

## SST Data File

In addition to the analytic SST modes the user can also specify an idealized aquaplanet SST pattern via the `_DOCN%AQPFILE_` option. The `aquapfile` namelist variable is used to specify the SST pattern in this mode. Note that this option has not been used or tested recently, so the user may experience difficulty trying to use this feature.
