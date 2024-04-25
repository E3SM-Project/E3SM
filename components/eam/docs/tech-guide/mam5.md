# Five-mode Modal Aerosol Model

## Overview

The Five-mode Modal Aerosol Model (MAM5) supersedes the MAM4 utilized in previous iterations of E3SM (E3SM-V1 and -V2). MAM5 introduces a fifth mode, specifically designed to represent stratospheric coarse mode aerosols, primarily originating from volcanic eruptions and DMS decomposition. This mode exclusively comprises sulfate aerosols, characterized by a smaller standard deviation (STD) value of 1.2. The STD value denotes the width of the aerosol mode, where a higher STD implies a greater gravitational settling effect (Wang et al., 2020; [@wang_aerosols_2020] Liu et al., 2012 [@liu_toward_2012]). By setting the STD to 1.2, the simulated properties of volcanic aerosols align closely with observational findings (Mills et al., 2016). [@mills_global_2016] MAM5 represents a pioneering aerosol model, effectively segregating tropospheric and stratospheric aerosols (Ke et al., in preparation), thereby mitigating the risk of overestimating dust and sea salt aerosols within the stratosphere in previous MAM4 (Visioni et al., 2021 [@visioni_limitations_2022]). Volcanic eruptions derived from Neely and Schmidt (2016) [@neely_iii_volcaneesm_2016].

## Namelist parameters

| Parameter                       | Description                                                                         | Default value               |
| ------------------------------- | ----------------------------------------------------------------------------------- | --------------------------- |
| `is_output_interactive_volc`    | Switch for diagnostic output of the stratospheric aerosol optics | `.false.`  |
