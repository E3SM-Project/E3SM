# Four-mode Modal Aerosol Module

## Overview

The representation of atmospheric aerosols and their roles in the Earth system by EAMv1/v2/v3 was inherited from the global aerosol-climate model EAMv0 and its four-mode modal aerosol module (MAM4), including Aitken, primary-carbon, accumulation, and coarse modes (Liu et al., 2016). [@liu_description_2016] It treats a combination of processes, controlling the evolution of aerosols that are either directly emitted or converted from precursor gases from a variety of natural and anthropogenic sources. The processes include transport (by grid-scale wind, subgrid turbulence, convection, and sedimentation), aerosol microphysics (i.e., particle nucleation, condensation/evaporation of trace gases, aging, and coagulation), cloud processing (i.e., aqueous chemistry, scavenging by hydrometeors, resuspension from evaporating hydrometeors, and wet deposition), and dry deposition. Aerosol species in the original MAM4 (Liu et al., 2016) [@liu_description_2016] include sulfate, primary organic aerosol (POA) or particulate organic matter (POM), secondary organic aerosol (SOA), black carbon (BC), sea salt, and mineral dust. As described by Wang et al. (2020), [@wang_aerosols_2020] the enhanced MAM4 in EAMv1/v2 added marine organic aerosol (MOA) to all four modes (Burrows et al., 2022). [@burrows_oceanfilms_2022] In MAM4 of EAMv3, the Aitken mode has sulfate, sea salt, SOA and MOA; the primary-carbon mode has BC, POA and MOA; the accumulation and coarse modes include all seven species. Ammonium (NH4) and nitrate (NO3) aerosols are also explicitly treated in EAMv3 (Wu et al., 2022), [@wu_development_2022] as an optional feature for research, in which new species (NH4, NO3, Ca, CO3, Na, Cl) are introduced to the Aitken, accumulation and coarse modes . All aerosol species within each of the four individual modes the MAM4 is assumed to be internally mixed and represented by a single number concentration, while particles are externally mixed among the different modes.

### Sea salt

In MAM4, sea salt aerosol is represented in the Aitken, accumulation, and coarse mode with particle emission size (diameter) ranges of 0.02-0.08, 0.08-1.0, and 1.0-10.0 μm, respectively. The emission flux of natural sea salt is first divided into 31 size bins, following the parameterization of Mårtensson et al. (2003) [@martensson_laboratory_2003] and Monahan et al. (1986), [@monahan_model_1986] and then redistributed to the three MAM4 size modes.

## Namelist parameters

| Parameter                | Description                                                                         | Default value               |
| ------------------------ | ----------------------------------------------------------------------------------- | --------------------------- |
| `mam_amicphys_optaa`     | Recommended option of the new time-splitting treatment of H2SO4 production and loss | `1` <!-- markdownlint-disable MD033 --><br> (0 to turn it off) |
| `n_so4_monolayers_pcage` | Number of monolayers required to age primary-carbon mode particles                  | `3`                         |
| `seasalt_emis_scale`     | Tuning parameter for sea salt emission                                              | `0.55`                      |
