# Five-mode Modal Aerosol Model

## Overview

The Five-mode Modal Aerosol Model (MAM5) supersedes the MAM4 utilized in previous iterations of E3SM (E3SM-V1 and -V2). MAM5 introduces a fifth mode, specifically designed to represent stratospheric coarse mode aerosols, primarily originating from volcanic eruptions and DMS decomposition. This mode exclusively comprises sulfate aerosols, characterized by a smaller standard deviation (STD) value of 1.2. The STD value denotes the width of the aerosol mode, where a higher STD implies a greater gravitational settling effect (Wang et al., 2019; Liu et al., 2012). By setting the STD to 1.2, the simulated properties of volcanic aerosols align closely with observational findings (Mills et al., 2016). MAM5 represents a pioneering aerosol model, effectively segregating tropospheric and stratospheric aerosols (Ke et al., in preparation), thereby mitigating the risk of overestimating dust and sea salt aerosols within the stratosphere in previous MAM4 (Visioni et al., 2021). Volcanic eruptions derived from Neely and Schmidt (2016).

Reference:

Liu, X., Easter, R. C., Ghan, S. J., Zaveri, R., Rasch, P., Shi, X., ... & Mitchell, D. (2012). Toward a minimal representation of aerosols in climate models: Description and evaluation in the Community Atmosphere Model CAM5. Geoscientific Model Development, 5(3), 709-739.

Mills, M. J., Schmidt, A., Easter, R., Solomon, S., Kinnison, D. E., Ghan, S. J., ... & Gettelman, A. (2016). Global volcanic aerosol properties derived from emissions, 1990–2014, using CESM1 (WACCM). Journal of Geophysical Research: Atmospheres, 121(5), 2332-2348.

Neely III, R. R., & Schmidt, A. (2016). VolcanEESM: Global volcanic sulphur dioxide (SO2) emissions database from 1850 to present.

Visioni, D., Tilmes, S., Bardeen, C., Mills, M., MacMartin, D. G., Kravitz, B., & Richter, J. (2021). Potential limitations of using a modal aerosol approach for sulfate geoengineering applications in climate models. Atmospheric Chemistry & Physics Discussions.

Wang, H., Easter, R. C., Zhang, R., Ma, P. L., Singh, B., Zhang, K., ... & Yoon, J. H. (2020). Aerosols in the E3SM Version 1: New developments and their impacts on radiative forcing. Journal of Advances in Modeling Earth Systems, 12(1), e2019MS001851.

## Namelist parameters
is_output_interactive_volc = .true.
output stratospheric AOD flag, default value = .false. 
