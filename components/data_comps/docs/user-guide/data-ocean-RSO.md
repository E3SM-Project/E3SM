# Data Ocean - Relaxed Slab Ocean (RSO)

The relaxed slab ocean (RSO) is similar in many ways to the [traditional slab ocean model](data-ocean-SOM.md), but uses a specified relaxation time scale to avoid the need for specified "Q-flux" data to represent the effects of ocean transport. The RSO implementation in E3SM was inspired by Zarzycki (2016)[@Zarzycki_TC-ocn-cpl_2016].

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
