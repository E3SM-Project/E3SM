# Rapid Radiative Transfer Model for GCMs

## Overview

The calculation of radiative energy flux through the atmosphere is done using the RRTMG radiation package (Iacono et al., 2008;[@iacono_radiative_2008] Mlawer et al., 1997[@mlawer_radiative_1997]). The details are consistent with the implementation in CAM5 described in Neale et al. (2012).[@neale_description_2012] Radiative fluxes are broadly split into shortwave and longwave and computed by separate codes. The shortwave solver uses the 2-stream approximation, while the longwave is an absorption/emission code. Both shortwave and longwave use the correlated k-distribution method for integration of fluxes. Subgrid cloud overlap is accounted for using the Monte Carlo Independent Column Approximation (MCICA; Pincus and Morcrette, 2003),[@pincus_fast_2003] assuming the cloudy portions of the column are maximally overlapped in vertically contiguous layers and randomly overlapped when two layers are separated by a completely clear layer. Cloud optics are parameterized as described in Neale et al.(2010).[@neale_description_2012]

## Namelist parameters

[RRTMG Namelist Parameters](../user-guide/namelist_parameters.md#rapid-radiative-transfer-model-for-gcms)
