# 1. Plant respiration

## 1.1 Autotrophic Respiration

The model treats maintenance and growth respiration processes separately (Lavigne and Ryan, 1997[@lavigne1997growth]; Sprugel, et al., 1995[@sprugel1995respiration]).

### 1.1.1 Maintenance Respiration

Maintenance respiration (MR) is a function of temperature and tissue N concentration (Ryan, 1991[@ryan1991effects]) for live biomass (excludes dead stem and coarse root pools) (Thornton and Rosenbloom, 2005[@thornton2005ecosystem]).
Rates for aboveground pools are based on the prognostic 2-meter air temperature, and rates for belowground pools (fine and coarse roots) depend on fractional rooting distribution with depth and the associated prognostic soil temperatures.

Implemented in subroutine: `MaintenanceResp`

$$
\begin{equation}
F^C_{mr,leaf} = (lmr_{sun} \times lai_{sun} + lmr_{sha} \times lai_{sha}) \times C_{umolCO2\_to\_gC}
\end{equation}
$$

$$
\begin{eqnarray}
\label{eqn_mr}
F^C_{mr,livestem} &=& S^N_{livestem} \times br_{mr} \times q_{10}^{\left(\frac{T_{2m} - T_{frz} - 20}{10}\right)} \\
F^C_{mr,livecroot} &=& S^N_{livecroot} \times br_{mr} \times q_{10}^{\left(\frac{T_{2m} - T_{frz} - 20}{10}\right)} \\
F^C_{mr,grain} &=& S^N_{grain} \times br_{mr} \times q_{10}^{\left(\frac{T_{2m} - T_{frz} - 20}{10}\right)} \\
F^C_{mr,fineroot} &=& \sum_{j=1} ^{nlevsoi} S^N_{fineroot} \times root_{fr,j} \times br_{mr} \times q_{10}^{\left(\frac{T_{s,j} - T_{frz} - 20}{10}\right)}
\end{eqnarray}
$$

$$
\begin{equation}
mr = F^C_{mr,leaf} + F^C_{mr,livestem} + F^C_{mr,livecroot} + F^C_{mr,grain} + F^C_{mr,fineroot}
\end{equation}
$$

### 1.1.2 Growth Respiration

Growth respiration is calculated as the parameter $gr_{perc}$, times the total carbon in new growth on a given timestep, based on construction costs for a range of woody and non-woody tissues (Larcher, 1995[@larcher2003physiological]).

Implemented in subroutine: `GrowthResp`

$$
\begin{eqnarray}
\label{eqn_gr}
F^C_{gr,leaf} &=& F^C_{alloc,leaf} \times gr_{perc} \\
F^C_{gr,leaf\_stor} &=& F^C_{alloc,leaf\_stor} \times gr_{perc} \times gr_{pnow} \\
F^C_{gr,froot} &=& F^C_{alloc,froot} \times gr_{perc} \\
F^C_{gr,froot\_stor} &=& F^C_{alloc,froot\_stor} \times gr_{perc} \times gr_{pnow} \\
F^C_{gr,livestem} &=& F^C_{alloc,livestem} \times gr_{perc} \\
F^C_{gr,livestem\_stor} &=& F^C_{alloc,livestem\_stor} \times gr_{perc} \times gr_{pnow} \\
F^C_{gr,deadstem} &=& F^C_{alloc,deadstem} \times gr_{perc} \\
F^C_{gr,deadstem\_stor} &=& F^C_{alloc,deadstem\_stor} \times gr_{perc} \times gr_{pnow} \\
F^C_{gr,livecroot} &=& F^C_{alloc,livecroot} \times gr_{perc} \\
F^C_{gr,livecroot\_stor} &=& F^C_{alloc,livecroot\_stor} \times gr_{perc} \times gr_{pnow} \\
F^C_{gr,deadcroot} &=& F^C_{alloc,deadcroot} \times gr_{perc} \\
F^C_{gr,deadcroot\_stor} &=& F^C_{alloc,deadcroot\_stor} \times gr_{perc} \times gr_{pnow} \\
F^C_{gr,grain} &=& F^C_{alloc,grain} \times gr_{perc} \\
F^C_{gr,grain\_stor} &=& F^C_{alloc,grain\_stor} \times gr_{perc} \times gr_{pnow} \\
\end{eqnarray}
$$

Parameter $gr_{pnow}$ is currently set to 1.0. This parameter could be changed to a smaller value to transfer some portion (1 - $gr_{pnow}$) of the growth respiration forward in time to occur at the time of growth display from storage.

## 1.2 Variable definitions

| Symbol  | Description |  Units| ELM Variable  |
| ----------- | ----------- | ----------- | ----------- |
| $lmr_{sun}$    | sunlit leaf maintenance respiration rate |  $umolCO_{2}~m^{-2}~sec^{-1}$ | `photosyns_vars%lmrsun_patch` |
| $lmr_{sha}$    | shaded leaf maintenance respiration rate |  $umolCO_{2}~m^{-2}~sec^{-1}$ | `photosyns_vars%lmrsha_patch` |
| $lai_{sun}$    | sunlit projected leaf area index |  - | `canopystate_vars%laisun_patch` |
| $lai_{sha}$    | shaded projected leaf area index |  - | `canopystate_vars%laisha_patch` |
| $br_{mr}$    | base rate for maintenance respiration |  $gC~gN^{-1}~sec^{-1}$ | `br_mr` |
| $q_{10}$   | temperature sensitivity for maintenance respiration |  - | `q10`  |
| $T_{2m}$   | 2 m height surface air temperature |  K | `veg_es%t_ref2m`  |
| $T_{s,j}$   | soil temperature at soil level j|  K | `col_es%t_soisno`  |
| $T_{frz}$   | freezing point of water (273.15 K) |  K | `SHR_CONST_TKFRZ`  |
| $root_{fr,j}$   |  fraction of roots in each soil layer |  - | `soilstate_vars%rootfr_patch`  |
| $S^N_{livestem}$    |  live stem N | $gN~m^{-2}$ | `veg_ns%livestemn`  |
| $S^N_{livecroot}$    |  live coarse root N | $gN~m^{-2}$ | `veg_ns%livecrootn`  |
| $S^N_{fineroot}$    |  fine root N | $gN~m^{-2}$ | `veg_ns%finerootn`  |
| $S^N_{grain}$    |  grain N | $gN~m^{-2}$ | `veg_ns%grainn`  |
| $F^C_{mr,leaf}$    |  leaf maintenance respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%leaf_mr`  |
| $F^C_{mr,livestem}$    |  live stem maintenance respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%livestem_mr`  |
| $F^C_{mr,livecroot}$    |  live coarse root maintenance respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%livecroot_mr`  |
| $F^C_{mr,fineroot}$    |  fine root maintenance respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%fineroot_mr`  |
| $F^C_{mr,grain}$    |  crop grain or organs maintenance respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%grain_mr`  |
| $gr_{perc}$    |  growth respiration parameter | - | `grperc`  |
| $gr_{pnow}$    |  growth respiration parameter | - | `grpnow`  |
| $F^C_{alloc,leaf}$ | allocation to leaf C | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_leafc` |
| $F^C_{alloc,leaf\_stor}$ | allocation to leaf C storage |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_leafc_storage` |
| $F^C_{alloc,froot}$ | allocation to fine root C | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_frootc` |
| $F^C_{alloc,froot\_stor}$ | allocation to fine root C storage | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_frootc_storage` |
| $F^C_{alloc,livestem}$ | allocation to live stem C | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_livestemc` |
| $F^C_{alloc,livestem\_stor}$ |  allocation to live stem C storage |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_livestemc_storage` |
| $F^C_{alloc,deadstem}$ | allocation to dead stem C |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_deadstemc` |
| $F^C_{alloc,deadstem\_stor}$ | allocation to dead stem C storage |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_deadstemc_storage` |
| $F^C_{alloc,livecroot}$ | allocation to live coarse root C |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_livecrootc` |
| $F^C_{alloc,livecroot\_stor}$ | allocation to live coarse root C storage |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_livecrootc_storage` |
| $F^C_{alloc,deadcroot}$ | allocation to dead coarse root C | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_deadcrootc` |
| $F^C_{alloc,deadcroot\_stor}$ | allocation to dead coarse root C storage |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_deadcrootc_storage` |
| $F^C_{alloc,grain}$ | allocation to grain C for prognostic crop | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_grainc` |
| $F^C_{alloc,grain\_stor}$ | allocation to grain C storage for prognostic crop |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_to_grainc_storage` |
| $F^C_{gr,leaf}$    |  leaf growth respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_leaf_gr`  |
| $F^C_{gr,leaf\_stor}$    |  leaf growth respiration to storage | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_leaf_storage_gr`  |
| $F^C_{gr,froot}$    |  fine root growth respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_froot_gr`  |
| $F^C_{gr,froot\_stor}$    |  fine root growth respiration to storage | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_froot_storage_gr`  |
| $F^C_{gr,livestem}$    |  livestem growth respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_livestem_gr`  |
| $F^C_{gr,livestem\_stor}$    |  livestem growth respiration to storage | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_livestem_storage_gr`  |
| $F^C_{gr,deadstem}$    |  deadstem growth respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_deadstem_gr`  |
| $F^C_{gr,deadstem\_stor}$    |  deadstem growth respiration to storage | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_deadstem_storage_gr`  |
| $F^C_{gr,livecroot}$    |  livecroot growth respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_livecroot_gr`  |
| $F^C_{gr,livecroot\_stor}$    |  livecroot growth respiration to storage | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_livecroot_storage_gr`  |
| $F^C_{gr,deadcroot}$    |  deadcroot growth respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_deadcroot_gr`  |
| $F^C_{gr,deadcroot\_stor}$    |  deadcroot growth respiration to storage | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_deadcroot_storage_gr`  |
| $F^C_{gr,grain}$    |  grain growth respiration | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_grain_gr`  |
| $F^C_{gr,grain\_stor}$    |  grain growth respiration to storage | $gC~m^{-2}~sec^{-1}$ | `veg_cf%cpool_grain_storage_gr`  |
