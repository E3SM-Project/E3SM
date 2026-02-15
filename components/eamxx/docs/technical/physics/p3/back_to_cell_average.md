# P3 Back To Cell Average

This document describes the `back_to_cell_average` bookkeeping step in P3 and
its unit-test coverage.

## Physical Formulation

P3 process rates are often computed for in-cloud or overlap regions. The model
state update uses cell-mean tendencies, so each tendency is mapped by the
appropriate cloud fraction:

$$
\left(\frac{\partial q}{\partial t}\right)_{\text{cell-avg}}
=
\left(\frac{\partial q}{\partial t}\right)_{\text{in-cloud}}
\times f_{\text{region}}
$$

Cloud overlap factors used in the implementation are:

$$
f_{ir} = \min(f_i, f_r),\quad
f_{il} = \min(f_i, f_l),\quad
f_{lr} = \min(f_l, f_r),\quad
f_{\text{glaciated}} = \max(10^{-4}, f_i - f_{il}).
$$

`f_glaciated` is used for the ice-only branch when
`use_separate_ice_liq_frac=true`.

## Implementation Details

- Source: `components/eamxx/src/physics/p3/impl/p3_back_to_cell_average_impl.hpp`.
- Runtime option: `use_separate_ice_liq_frac`.
- Threshold: `cld_frac_glaciated` has a hard minimum of `1e-4`.
- Two nucleation tendencies are currently passed through unchanged:
  - `qv2qi_nucleat_tend`
  - `ni_nucleat_tend`

## Process Mapping

| Process | Variable | Cloud Fraction | Rationale |
| --- | --- | --- | --- |
| Accretion cloud->rain mass | `qc2qr_accret_tend` | $f_{lr}$ | Requires liquid-rain overlap |
| Rain evaporation mass | `qr2qv_evap_tend` | $f_r$ | Occurs in rain region |
| Autoconversion mass | `qc2qr_autoconv_tend` | $f_l$ | Occurs in liquid cloud |
| Accretion cloud number | `nc_accret_tend` | $f_{lr}$ | Number sink tied to accretion overlap |
| Cloud self-collection number | `nc_selfcollect_tend` | $f_l$ | Liquid-cloud process |
| Autoconversion cloud number transfer | `nc2nr_autoconv_tend` | $f_l$ | Liquid-cloud process |
| Rain self-collection number | `nr_selfcollect_tend` | $f_r$ | Rain-cloud process |
| Rain evaporation number sink | `nr_evap_tend` | $f_r$ | Rain-cloud process |
| Autoconversion rain number source | `ncautr` | $f_{lr}$ | Liquid-rain overlap |
| Ice sublimation mass | `qi2qv_sublim_tend` | $f_{\text{glaciated}}$ or $f_i$ | Runtime branch |
| Ice shedding rain number source | `nr_ice_shed_tend` | $f_{il}$ | Mixed-phase interaction |
| Heterogeneous freezing mass | `qc2qi_hetero_freeze_tend` | $f_{il}$ | Mixed-phase interaction |
| Rain collected by ice mass | `qr2qi_collect_tend` | $f_{ir}$ | Ice-rain overlap |
| Ice shedding rain mass source | `qc2qr_ice_shed_tend` | $f_{il}$ | Mixed-phase interaction |
| Ice melting mass | `qi2qr_melt_tend` | $f_i$ | Ice-cloud process |
| Cloud collected by ice mass | `qc2qi_collect_tend` | $f_{il}$ | Mixed-phase interaction |
| Rain immersion freezing mass | `qr2qi_immers_freeze_tend` | $f_r$ | Rain-cloud process |
| Ice melting rain number source | `ni2nr_melt_tend` | $f_i$ | Ice-cloud process |
| Cloud collected by ice number | `nc_collect_tend` | $f_{il}$ | Mixed-phase interaction |
| Shedding cloud number sink | `ncshdc` | $f_{il}$ | Mixed-phase interaction |
| Cloud immersion-freezing number | `nc2ni_immers_freeze_tend` | $f_l$ | Implemented as liquid-cloud process |
| Rain collected by ice number | `nr_collect_tend` | $f_{ir}$ | Ice-rain overlap |
| Ice self-collection number | `ni_selfcollect_tend` | $f_i$ | Ice-cloud process |
| Vapor deposition to ice mass | `qv2qi_vapdep_tend` | $f_{\text{glaciated}}$ or $f_i$ | Runtime branch |
| Rain immersion-freezing number | `nr2ni_immers_freeze_tend` | $f_r$ | Rain-cloud process |
| Ice sublimation number sink | `ni_sublim_tend` | $f_{\text{glaciated}}$ or $f_i$ | Runtime branch |
| Ice nucleation mass | `qv2qi_nucleat_tend` | unchanged | Already cell-averaged |
| Ice nucleation number | `ni_nucleat_tend` | unchanged | Already cell-averaged |
| Bergeron mass transfer | `qc2qi_berg_tend` | $f_{il}$ | Mixed-phase process |
| Het-freezing number counter | `ncheti_cnt` | $f_l$ | Counter mapped over liquid cloud |
| Het-freezing mass counter | `qcheti_cnt` | $f_l$ | Counter mapped over liquid cloud |
| Contact nucleation number counter | `nicnt` | $f_l$ | Counter mapped over liquid cloud |
| Contact nucleation mass counter | `qicnt` | $f_l$ | Counter mapped over liquid cloud |
| Nucleation number counter | `ninuc_cnt` | $f_l$ | Counter mapped over liquid cloud |
| Nucleation mass counter | `qinuc_cnt` | $f_l$ | Counter mapped over liquid cloud |

## Property Tests

The unit tests live in
`components/eamxx/src/physics/p3/tests/p3_back_to_cell_average_unit_tests.cpp`.

### Test Organization

The Catch2 test case `p3_back_to_cell_average` is split into independent
sections:

- `process_location`
- `bounds_and_scaling`
- `cloud_fraction_logic`
- `runtime_option`
- `edge_cases`
- `bookkeeping_diagnostics` (non-gating)
- `bfb`

### Tolerance Philosophy

| Test Category | Tolerance Type | Value | Rationale |
| --- | --- | --- | --- |
| Process-location mapping | Identity | `10 * epsilon` | Exact implementation formula |
| Runtime branch mapping | Identity | `10 * epsilon` | Exact branch behavior |
| Intersection consistency | Identity | `10 * epsilon` | `min`/`max` identity checks |
| Scaling reduction | Identity | `10 * epsilon` | For factors in `[0,1]` |
| Non-negativity | Absolute | `1e-30` | Numerical zero proxy |
| Glaciated threshold | Hard bound | `>= 1e-4` | Explicit implementation floor |

### Parameter Sweep Strategy

A 3D linear sweep over cloud fractions is used:

- $f_i, f_l, f_r \in [0,1]$
- `n_frac = 15` per dimension
- `num_cases = 15^3 = 3375`

All tendency inputs are deterministic, distinct, and non-negative.

### Number-Tendency Validation Policy

The suite does **not** enforce a global number-conservation identity
($\sum \partial N / \partial t = 0$). In P3, self-collection, evaporation/
sublimation, nucleation, and shedding are not globally number-conservative.

Instead, number tendencies are validated via process-location mapping checks:
for each number tendency, the test verifies the exact cloud-fraction factor used
by the implementation.

### Bookkeeping Diagnostics

A non-gating diagnostic reports the maximum residual of an assembled total-water
bookkeeping sum over the sweep. This is informational and not used as a pass/
fail criterion.
