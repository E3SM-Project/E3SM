# 1. Carbon, nitrogen, and phosphorus allocation

E3SM is fully prognostic for all carbon, nitrogen, and phosphorus state variables in the vegetation, litter, and soil organic matter pools.
At each model time step and for each PFT occupying a soil column, the carbon available for allocation to new growth ($F^C_{avail}$) is calculated. The carbon available is then allocated to new plant growth.

## 1.1 Carbon available for new growth

Implemented in subroutine: `Allocation1_PlantNPDemand`

$$
\begin{equation}
\label{eqn_gpp}
F^C_{GPP} = A^{sun} + A^{shade}
\end{equation}
$$

$A^{sun}$ - carbon fixation from sunlight canopy
$A^{shade}$ - carbon fixation from shaded canopy

On each model time step and for each PFT sharing space on a soil column,
the carbon available for allocation to new growth ($F^C_{avail}$) is calculated as:

$$
\begin{equation}
F^C_{avail} = max\left( 0, F^C_{GPP} - F^C_{MR} \right) \\
\end{equation}
$$

$$
\begin{equation}
\label{eqn_Fc_storage}
F^C_{storage\_recov} =
\begin{cases}
-\frac{F^C_{storage}}{86400\tau _{xs}} \qquad \text{for} \qquad F^C_{storage} < 0 \\
F^C_{avail} \qquad \text{for} \qquad \left(F^C_{storage\_recov} > F^C_{avail} \right) \\
\end{cases}
\end{equation}
$$

$$
\begin{equation}
\label{eqn_Fc_avail}
F^C_{avail} =
\begin{cases}
F^C_{avail} - F^C_{storage\_recov} \qquad \text{for} \qquad \left(F^C_{storage\_recov} < F^C_{avail} \right)  \\
0 \qquad \text{for} \qquad \left(F^C_{storage\_recov} > F^C_{avail} \right)  \\
\end{cases}
\end{equation}
$$

## 1.2 Plant nitrogen and phosphorus demand for new growth

Implemented in subroutine: `Allocation1_PlantNPDemand`

Total plant nitrogen and phosphorus demand for each time step is calculated from carbon allometry and nitrogen and phosphorus concentrations for each tissue type (specified by PFT).

$$
\begin{equation}
\label{eqn_Ndemand}
F^N_{plant\_demand} = F^C_{avail} \frac{N_{allom}}{C_{allom}}
\end{equation}
$$

$$
\begin{equation}
\label{eqn_Pdemand}
F^P_{plant\_demand} = F^C_{avail} \frac{P_{allom}}{C_{allom}}
\end{equation}
$$

$$
\begin{equation}
\label{eqn_Callom}
C_{allom} = \begin{cases}
{\left(1+g_{1} \right)\left(1+a_{1} +a_{3} \left(1+a_{2} \right)\right)\qquad \text{for woody PFT}} \\
{1+g_{1} +a_{1} \left(1+g_{1} \right)\qquad \qquad \text{for non-woody PFT}} \\
{\left(1+g_{1} \right)\left(1+a_{1}+a_{5} +a_{3} \left(1+a_{2} \right)\right)\qquad \text{for crop PFT}}
\end{cases}
\end{equation}
$$

$$
\begin{equation}
\label{eqn_Nallom}
N_{allom} = \begin{cases}
{\frac{1}{CN_{leaf} } +\frac{a_{1} }{CN_{fr} } +\frac{a_{3} a_{4} \left(1+a_{2} \right)}{CN_{lw} } + \frac{a_{3} \left(1-a_{4} \right)\left(1+a_{2} \right)}{CN_{dw} } \qquad \text{for woody PFT}} \\
{\frac{1}{CN_{leaf} } +\frac{a_{1} }{CN_{fr} } \qquad \qquad \qquad \qquad \qquad \text{for non-woody PFT}} \\
{\frac{1}{CN_{leaf} } +\frac{a_{1} }{CN_{fr} } +\frac{a_{5} }{CN_g} +\frac{a_{3} a_{4} \left(1+a_{2} \right)}{CN_{lw} } + \frac{a_{3} \left(1-a_{4} \right)\left(1+a_{2} \right)}{CN_{dw} } \qquad \text{for crop PFT}}
\end{cases}
\end{equation}
$$

$$
\begin{equation}
\label{eqn_Pallom}
P_{allom} =\begin{cases}
{\frac{1}{CP_{leaf} } +\frac{a_{1} }{CP_{fr} } +\frac{a_{3} a_{4} \left(1+a_{2} \right)}{CP_{lw} } + \frac{a_{3} \left(1-a_{4} \right)\left(1+a_{2} \right)}{CP_{dw} } \qquad \text{for woody PFT}} \\
{\frac{1}{CP_{leaf} } +\frac{a_{1} }{CP_{fr} } \qquad \qquad \qquad \qquad \qquad \text{for non-woody PFT}} \\
{\frac{1}{CP_{leaf} } +\frac{a_{1} }{CP_{fr} } +\frac{a_{5} }{CP_g} +\frac{a_{3} a_{4} \left(1+a_{2} \right)}{CP_{lw} } + \frac{a_{3} \left(1-a_{4} \right)\left(1+a_{2} \right)}{CP_{dw} } \qquad \text{for crop PFT}}
\end{cases}
\end{equation}
$$

The allometric parameters relate allocation between various tissue types:

| Coefficient | Description |
| ----------- | ----------- |
| $a_{1}$   |  ratio of new fine root : new leaf carbon allocation |
| $a_{2}$   |  ratio of new coarse root : new stem carbon allocation |
| $a_{3}$   |  ratio of new stem : new leaf carbon allocation |
| $a_{4}$   |  ratio of new live wood : new total wood carbon allocation |
| $a_{5}$   |  ratio of grain : new total carbon allocation |
| $g_{1}$   |  ratio of growth respiration carbon : new growth carbon |

Parameters $a_{1}$, $a_{2}$, $a_{4}$, and $a_{5}$ are defined as constants for various non-crop PFTs (Table ) and estimated for crop pfts.
$g_{1}$ = 0.3 (unitless) for non-crop PFTS and 0.25 or 0.11 for crop PFTS.

A dynamic allocation scheme is used for woody vegetation.

$$
\begin{equation}
a_{3} = max \left( \left( \frac{2.7}{1+e^{-0.004NPP_{ann} -300} } -0.4 \right), 0.2  \right)
\end{equation}
$$

where $NPP_{ann}$ is the annual sum of NPP from the previous year.

## 1.3 Carbon, nitrogen, and phosphorus allocation to new growth

Implemented in subroutine: `Allocation3_PlantCNPAlloc`

### 1.3.1 Competition scaled by relative nutrient demand (RD) (CTC - Converging Trophic Cascade)

based on Thornton et al. 2002[@thornton2002modeling] and Yang et al., 2014[@yang2014role]

$$
\begin{eqnarray}
\label{eqn_NP_alloc}
    F^N_{alloc} &=& F^N_{plant\_demand} + N_{retrans} \\
    F^P_{alloc} &=& F^P_{plant\_demand} + P_{retrans}
\end{eqnarray}
$$

$$
\label{eqn_CNP_alloc}
\begin{array}{l}
if \left( F^N_{alloc} \frac{C_{allom}}{N_{allom}} < F^P_{alloc} \frac{C_{allom}}{P_{allom}} \right) then
    \begin{cases}
        F^C_{alloc} = F^N_{alloc} \frac{C_{allom}}{N_{allom}} \\
        F^P_{alloc} = F^N_{alloc} \frac{P_{allom}}{N_{allom}} \\
    \end{cases} \\
if \left( F^N_{alloc} \frac{C_{allom}}{N_{allom}} > F^P_{alloc} \frac{C_{allom}}{P_{allom}} \right) then
    \begin{cases}
        F^C_{alloc} = F^P_{alloc} \frac{C_{allom}}{P_{allom}} \\
        F^N_{alloc} = F^P_{alloc} \frac{N_{allom}}{P_{allom}} \\
    \end{cases} \\
\end{array}
$$

$$
\begin{eqnarray}
\label{eqn_CNP_alloc_final}
F^C_{alloc} &=& \frac{F^C_{alloc}}{C_{allom}} \\
F^N_{alloc} &=& \frac{F^N_{alloc}}{N_{allom}} \\
F^P_{alloc} &=& \frac{F^P_{alloc}}{P_{allom}} \\
\end{eqnarray}
$$

The allocation fluxes of carbon to display and storage pools for the various tissue types are given as:

$$
\begin{eqnarray}
\label{eqn_C_newgrowth}
F^C_{alloc,leaf} &=& F^C_{alloc} f_{cur} \\
F^C_{alloc,leaf\_stor} &=& F^C_{alloc} \left(1-f_{cur} \right) \\
F^C_{alloc,froot} &=& F^C_{alloc} a_{1} f_{cur} \\
F^C_{alloc,froot\_stor} &=& F^C_{alloc} a_{1} \left(1-f_{cur} \right) \\
F^C_{alloc,livestem} &=& F^C_{alloc} a_{3} a_{4} f_{cur} \\
F^C_{alloc,livestem\_stor} &=& F^C_{alloc} a_{3} a_{4} \left(1-f_{cur} \right) \\
F^C_{alloc,deadstem} &=& F^C_{alloc} a_{3} \left(1-a_{4} \right)f_{cur} \\
F^C_{alloc,deadstem\_stor} &=& F^C_{alloc} a_{3} \left(1-a_{4} \right)\left(1-f_{cur} \right) \\
F^C_{alloc,livecroot} &=& F^C_{alloc} a_{2} a_{3} a_{4} f_{cur} \\
F^C_{alloc,livecroot\_stor} &=& F^C_{alloc} a_{2} a_{3} a_{4} \left(1-f_{cur} \right) \\
F^C_{alloc,deadcroot} &=& F^C_{alloc} a_{2} a_{3} \left(1-a_{4} \right)f_{cur} \\
F^C_{alloc,deadcroot\_stor} &=& F^C_{alloc} a_{2} a_{3} \left(1-a_{4} \right)\left(1-f_{cur} \right) \\
F^C_{alloc,grain} &=& F^C_{alloc} a_{5} f_{cur} \\
F^C_{alloc,grain\_stor} &=& F^C_{alloc} a_{5} \left(1-f_{cur} \right)
\end{eqnarray}
$$

The allocation fluxes of nitrogen to display and storage pools for the various tissue types are given as:

$$
\begin{eqnarray}
\label{eqn_N_newgrowth}
F^N_{demand,leaf} &=& \frac{F^N_{alloc} }{CN_{leaf} } f_{cur} \\
F^N_{demand,leaf\_stor} &=& \frac{F^N_{alloc} }{CN_{leaf} } \left(1-f_{cur} \right) \\
F^N_{demand,froot} &=& \frac{F^N_{alloc} a_{1} }{CN_{fr} } f_{cur} \\
F^N_{demand,froot\_stor} &=& \frac{F^N_{alloc} a_{1} }{CN_{fr} } \left(1-f_{cur} \right) \\
F^N_{demand,livestem} &=& \frac{F^N_{alloc} a_{3} a_{4} }{CN_{lw} } f_{cur} \\
F^N_{demand,livestem\_stor} &=& \frac{F^N_{alloc} a_{3} a_{4} }{CN_{lw} } \left(1-f_{cur} \right) \\
F^N_{demand,deadstem} &=& \frac{F^N_{alloc} a_{3} \left(1-a_{4} \right)}{CN_{dw} } f_{cur} \\
F^N_{demand,deadstem\_stor} &=& \frac{F^N_{alloc} a_{3} \left(1-a_{4} \right)}{CN_{dw} } \left(1-f_{cur} \right) \\
F^N_{demand,livecroot} &=& \frac{F^N_{alloc} a_{2} a_{3} a_{4} }{CN_{lw} } f_{cur} \\
F^N_{demand,livecroot\_stor} &=& \frac{F^N_{alloc} a_{2} a_{3} a_{4} }{CN_{lw} } \left(1-f_{cur} \right) \\
F^N_{demand,deadcroot} &=& \frac{F^N_{alloc} a_{2} a_{3} \left(1-a_{4} \right)}{CN_{dw} } f_{cur} \\
F^N_{demand,deadcroot\_stor} &=& \frac{F^N_{alloc,leaf} a_{2} a_{3} \left(1-a_{4} \right)}{CN_{dw} } \left(1-f_{cur} \right) \\
\end{eqnarray}
$$

$$
\begin{eqnarray}
F^N_{demand,grain} &=& \frac{F^N_{alloc} a_{5} }{CN_{g} } f_{cur} \\
F^N_{demand,grain\_stor} &=& \frac{F^N_{alloc} a_{5} }{CN_{g} } \left(1-f_{cur} \right)
\end{eqnarray}
$$

The allocation fluxes of phosphorus to display and storage pools for the various tissue types are given as:

$$
\begin{eqnarray}
\label{eqn_P_newgrowth}
F^P_{demand,leaf} &=& \frac{F^P_{alloc} }{CP_{leaf} } f_{cur} \\
F^P_{demand,leaf\_stor} &=& \frac{F^P_{alloc} }{CP_{leaf} } \left(1-f_{cur} \right) \\
F^P_{demand,froot} &=& \frac{F^P_{alloc} a_{1} }{CP_{fr} } f_{cur} \\
F^P_{demand,froot\_stor} &=& \frac{F^P_{alloc} a_{1} }{CP_{fr} } \left(1-f_{cur} \right) \\
F^P_{demand,livestem} &=& \frac{F^P_{alloc} a_{3} a_{4} }{CP_{lw} } f_{cur} \\
F^P_{demand,livestem\_stor} &=& \frac{F^P_{alloc} a_{3} a_{4} }{CP_{lw} } \left(1-f_{cur} \right) \\
F^P_{demand,deadstem} &=& \frac{F^P_{alloc} a_{3} \left(1-a_{4} \right)}{CP_{dw} } f_{cur} \\
F^P_{demand,deadstem\_stor} &=& \frac{F^P_{alloc} a_{3} \left(1-a_{4} \right)}{CP_{dw} } \left(1-f_{cur} \right) \\
F^P_{demand,livecroot} &=& \frac{F^P_{alloc} a_{2} a_{3} a_{4} }{CP_{lw} } f_{cur} \\
F^P_{demand,livecroot\_stor} &=& \frac{F^P_{alloc} a_{2} a_{3} a_{4} }{CP_{lw} } \left(1-f_{cur} \right) \\
F^P_{demand,deadcroot} &=& \frac{F^P_{alloc} a_{2} a_{3} \left(1-a_{4} \right)}{CP_{dw} } f_{cur} \\
F^P_{demand,deadcroot\_stor} &=& \frac{F^P_{alloc,leaf} a_{2} a_{3} \left(1-a_{4} \right)}{CP_{dw} } \left(1-f_{cur} \right) \\
\end{eqnarray}
$$

$$
\begin{eqnarray}
F^P_{demand,grain} &=& \frac{F^P_{alloc} a_{5} }{CP_{g} } f_{cur} \\
F^P_{demand,grain\_stor} &=& \frac{F^P_{alloc} a_{5} }{CP_{g} } \left(1-f_{cur} \right)
\end{eqnarray}
$$

### 1.3.2 ECA - competition scaled by Equilibrium Chemical Approximation (ECA) theory

Zhu et al., 2019 [@zhu2019representing]

### 1.3.3 MIC - always assume that soil microbes out-compete plants

## 1.4 Variable definitions

| Symbol  | Description |  Units| ELM Variable  |
| ----------- | ----------- | ----------- | ----------- |
| $A^{sun}$   |  C fixation from sunlit canopy |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%psnsun_to_cpool` |
| $A^{shade}$   | C fixation from shaded canopy |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%psnshade_to_cpool` |
| $F^C_{GPP}$   |  GPP flux before downregulation |  $gC~m^{-2}~sec^{-1}$ | `veg_cf%gpp_before_downreg` |
| $F^C_{MR}$   |  maintenance respiration | $gC~m^{-2}~sec^{-1}$ | mr |
| $F^C_{avail}$  | C flux available for allocation | $gC~m^{-2}~sec^{-1}$ | `veg_cf%availc` |
| $C_{allom}$ |  C allocation index | -  |`cnstate_vars%c_allometry_patch` |
| $N_{allom}$ | N allocation index | -  |`cnstate_vars%n_allometry_patch` |
| $P_{allom}$ | P allocation index | -  |`cnstate_vars%p_allometry_patch` |
| $F^N_{plant\_demand}$ | N flux required to support initial GPP | $gN~m^{-2}~sec^{-1}$ | `veg_nf%plant_ndemand` |
| $F^P_{plant\_demand}$ | P flux required to support initial GPP |  $gP~m^{-2}~sec^{-1}$ | `veg_pf%plant_pdemand` |
| $CN_{leaf}$ | leaf C:N | $gC~gN^{-1}$ | `veg_vp%leafcn` |
| $CN_{fr}$ | fine root C:N |  $gC~gN^{-1}$ | `veg_vp%frootcn` |
| $CN_{lw}$ | live wood C:N | $gC~gN^{-1}$ | `veg_vp%livewdcn` |
| $CN_{dw}$ | dead wood C:N | $gC~gN^{-1}$ | `veg_vp%deadwdcn` |
| $CN_{g}$ | grain C:N | $gC~gN^{-1}$ | `veg_vp%graincn` |
| $CP_{leaf}$ | leaf C:P | $gC~gP^{-1}$ | `veg_vp%leafcp` |
| $CP_{fr}$ | fine root C:P | $gC~gP^{-1}$ | `veg_vp%frootcp` |
| $CP_{lw}$ | live wood C:P | $gC~gP^{-1}$ | `veg_vp%livewdcp` |
| $CP_{dw}$ | dead wood C:P | $gC~gP^{-1}$ | `veg_vp%deadwdcp` |
| $CP_{g}$ | grain C:P |  $gC~gP^{-1}$ | `veg_vp%graincp` |
| $N_{retrans}$ | deployment of retranslocated N | $gN~m^{-2}~sec^{-1}$ | `veg_nf%retransn_to_npool` |
| $P_{retrans}$ | deployment of retranslocated P |  $gP~m^{-2}~sec^{-1}$ | `veg_pf%retransp_to_ppool` |
| $F^C_{alloc}$ | total allocated C flux | $gC~m^{-2}~sec^{-1}$ | `veg_cf%plant_calloc` |
| $F^N_{alloc}$ | total allocated N flux | $gN~m^{-2}~sec^{-1}$ | `veg_nf%plant_nalloc` |
| $F^P_{alloc}$ | total allocated P flux | $gP~m^{-2}~sec^{-1}$ | `veg_pf%plant_palloc` |
| $f_{cur}$ | fraction of allocation that goes to currently displayed growth, remainder to storage | - | `veg_vp%fcur` |
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
| $F^N_{alloc,leaf}$ | allocation to leaf N | $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_leafn` |
| $F^N_{alloc,leaf\_stor}$ | allocation to leaf N storage |  $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_leafn_storage` |
| $F^N_{alloc,froot}$ | allocation to fine root N | $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_frootn` |
| $F^N_{alloc,froot\_stor}$ | allocation to fine root N storage | $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_frootn_storage` |
| $F^N_{alloc,livestem}$ | allocation to live stem N | $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_livestemn` |
| $F^N_{alloc,livestem\_stor}$ |  allocation to live stem N storage |  $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_livestemn_storage` |
| $F^N_{alloc,deadstem}$ | allocation to dead stem N |  $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_deadstemn` |
| $F^N_{alloc,deadstem\_stor}$ | allocation to dead stem N storage |  $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_deadstemn_storage` |
| $F^N_{alloc,livecroot}$ | allocation to live coarse root N |  $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_livecrootn` |
| $F^N_{alloc,livecroot\_stor}$ | allocation to live coarse root N storage |  $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_livecrootn_storage` |
| $F^N_{alloc,deadcroot}$ | allocation to dead coarse root N | $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_deadcrootn` |
| $F^N_{alloc,deadcroot\_stor}$ | allocation to dead coarse root N storage |  $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_deadcrootn_storage` |
| $F^N_{alloc,grain}$ | allocation to grain N for prognostic crop | $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_grainn` |
| $F^N_{alloc,grain\_stor}$ | allocation to grain N storage for prognostic crop |  $gN~m^{-2}~sec^{-1}$ | `veg_nf%npool_to_grainn_storage` |
| $F^P_{alloc,leaf}$ | allocation to leaf P | $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_leafp` |
| $F^P_{alloc,leaf\_stor}$ | allocation to leaf P storage |  $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_leafp_storage` |
| $F^P_{alloc,froot}$ | allocation to fine root P | $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_frootp` |
| $F^P_{alloc,froot\_stor}$ | allocation to fine root P storage | $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_frootp_storage` |
| $F^P_{alloc,livestem}$ | allocation to live stem P | $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_livestemp` |
| $F^P_{alloc,livestem\_stor}$ |  allocation to live stem P storage |  $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_livestemp_storage` |
| $F^P_{alloc,deadstem}$ | allocation to dead stem P |  $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_deadstemp` |
| $F^P_{alloc,deadstem\_stor}$ | allocation to dead stem P storage |  $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_deadstemp_storage` |
| $F^P_{alloc,livecroot}$ | allocation to live coarse root P |  $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_livecrootp` |
| $F^P_{alloc,livecroot\_stor}$ | allocation to live coarse root P storage |  $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_livecrootp_storage` |
| $F^P_{alloc,deadcroot}$ | allocation to dead coarse root P | $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_deadcrootp` |
| $F^P_{alloc,deadcroot\_stor}$ | allocation to dead coarse root P storage |  $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_deadcrootp_storage` |
| $F^P_{alloc,grain}$ | allocation to grain P for prognostic crop | $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_grainp` |
| $F^P_{alloc,grain\_stor}$ | allocation to grain P storage for prognostic crop |  $gP~m^{-2}~sec^{-1}$ | `veg_pf%ppool_to_grainp_storage` |
