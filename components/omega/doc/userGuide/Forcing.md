(omega-user-forcing)=

# Forcing

This page documents the user-facing configuration and behavior for current forcing in Omega:

- Surface stress forcing (e.g. wind stress)
- Coupled tracer flux forcing (mass, energy and salt)
- Surface tracer restoring

## Surface stress forcing

Surface stress forcing adds momentum tendency from surface stress (e.g. wind).

(omega-user-forcing-sfc-stress)=
### Surface stress forcing configuration

Surface stress forcing behavior is controlled by two configuration blocks:

```yaml
Omega:
  SfcStress:
    InterpType: Isotropic

  Tendencies:
    SfcStressForcingTendencyEnable: true
```

- `SfcStress.InterpType`
  - `Isotropic`: isotropic cell-to-edge interpolation for surface stress
  - `Anisotropic`: anisotropic interpolation option
- `Tendencies.SfcStressForcingTendencyEnable`: switch to enable surface stress forcing tendency

See {ref}`Interpolate cells to edges <omega-user-horz-operators-interp>` for a detailed description of the cell-to-edge interpolation operators (Isotropic and Anisotropic) used by SfcStress.InterpType.

### Required input fields

Surface stress forcing uses surface stress input fields:

- `SfcStressZonal`
- `SfcStressMeridional`

These are stored in forcing variables and used to form edge-normal stress
(`NormalStressEdge`) that enters momentum tendencies.

## Surface thickness and tracer flux forcing

Surface thickness and tracer flux forcing applies ocean-atmosphere and ocean-sea ice fluxes from the other model
components (atmosphere, sea ice) to the thickness and tracer equations. This enables
the ocean to respond to heat, freshwater, and salt exchanges at the surface. These fluxes can be from data or (active) coupled components.

### Surface thickness and tracer flux forcing configuration

Surface thickness and tracer flux forcing is controlled by two configuration flags:

```yaml
Omega:
  Tendencies:
    SfcThicknessForcingTendencyEnable: false
    SfcTracerForcingTendencyEnable: false
```

- `Tendencies.SfcThicknessForcingTendencyEnable`: enables coupled freshwater and salt flux forcing on thickness
- `Tendencies.SfcTracerForcingTendencyEnable`: enables coupled heat and salt flux forcing on tracers


### Required input fields

Coupled tracer flux forcing uses 13 auxiliary fields organized by type:

**Freshwater mass fluxes (kg m⁻² s⁻¹):**
- `SnowFlux`: precipitation from snow
- `RainFlux`: precipitation from rain
- `EvaporationFlux`: evaporative water loss
- `SeaIceFreshWaterFlux`: freshwater mass flux from sea-ice melt or formation
- `IceRunoffFlux`: runoff from land ice
- `RiverRunoffFlux`: runoff from rivers

**Heat/Enthalpy fluxes (W m⁻²):**
- `LatentHeatFluxEvap`: latent heat flux from evaporation phase change
- `SensibleHeatFlux`: sensible heat flux
- `LongWaveHeatFluxUp`: upward longwave radiation
- `LongWaveHeatFluxDown`: downward longwave radiation
- `SeaIceHeatFlux`: heat/energy from sea-ice interaction (incl. enthalpy of meltwater)
- `ShortWaveHeatFlux`: shortwave (solar) radiation

**Salt mass flux (kg m⁻² s⁻¹):**
- `SeaIceSaltFlux`: salt flux from sea-ice formation/melt processes

These fields are populated by external coupling components (typically atmosphere
and ice models). Omega assumes the incoming values match the documented units.
For now, there are assumed to come from a `forcing.nc` file, but later will be provided
by the equivalent `ocn_comp_mct.F`.

### Notes

- Coupled fluxes are applied only at the surface layer (top active layer) for each cell. In the future, vertical spreading of contributions from river runoff will be needed.
- Pseudo-thickness tendency is computed from the (six) freshwater mass fluxes and the salt mass flux
  `SeaIceSaltFlux`, converted to a pseudo-thickness change.
- Temperature tendency is computed from direct heat flux plus
  mass-flux enthalpy terms, converted to conservative-temperature tendency via
  $H_{\text{FluxFac}} = 1.0 / (\rho_{sw} c^0_{p,sw})$ where $c^0_{p,sw}$ is the reference
  specific heat of seawater defined by TEOS-10. The enthalpy associated with mass fluxes is currently hard-coded to SST for liquid fluxes and the freezing temperature for solid fluxes (which are melted using a constant latent heat of fusion). Note that the enthalpy of liquid meltwater from sea ice is already included in `SeaIceHeatFlux`.
- Salinity tendency from `SeaIceSaltFlux` is scaled by
  $S_{\text{FluxFac}} = 1.0e3 / \rho_{sw}$ to account for unit conversion from
  kg/(m²·s) to salinity units (g/kg).
- Fluxes are assumed to be in the documented units (i.e. net mass fluxes);
  any unit conversion should be performed by the coupling component before providing flux
  values to Omega.
- The reference density used here ($\rho_{sw}$) is not a Boussinesq density, it is the
  conversion factor from mass to pseudo-thickness.
- No iceberg fluxes are included for now.

## Penetrating shortwave radiation

Penetrating shortwave radiation distributes incoming solar heat through the
active water column instead of applying all shortwave heating to the surface
layer. It is enabled with:

```yaml
Omega:
  Tendencies:
    SfcTracerForcingTendencyEnable: true
    PenetratingShortwaveTendency:
      Enable: true
      NearIrFraction: 0.58
      NearIrCoeff: 2.86
      RedFraction: 0.21
      BlueFraction: 0.21

  IOStreams:
    ShortwaveExtinctionIn:
      Filename: shortwave_extinction_coeffs.nc
```

The `ShortwaveExtinctionIn` stream is read at startup. Its file must provide
the cell-centered fields `ExtinctionCoeffRedCell` and
`ExtinctionCoeffBlueCell`, both in `m^-1`. These fields provide the red-band
($k_r$) and blue/green-band ($k_b$) attenuation coefficients for the annual
average used by the simulation.

The incoming `ShortWaveHeatFlux` is split into three bands:

- near infrared: fraction `NearIrFraction`, extinction coefficient
  `NearIrCoeff` ($0.58$ and $2.86\ \mathrm{m}^{-1}$ by default)
- red: fraction `RedFraction` ($0.21$ by default), extinction coefficient $k_r$
- blue/green: fraction `BlueFraction` ($0.21$ by default), extinction coefficient $k_b$

The fractions (`NearIrFraction`, `RedFraction`, and `BlueFraction`) and the
near-infrared extinction coefficient (`NearIrCoeff`) can be configured under
`Omega.Tendencies.PenetratingShortwaveTendency`. By default, the red and blue
fractions are set to 0.21 as in [Manizza et al. (2005)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2004GL020778).
The sum of `NearIrFraction`, `RedFraction`, and `BlueFraction` must equal 1.

The attenuated flux is converted to temperature heating from the difference
between the flux at the top and bottom of each layer. Any flux remaining at the
bottom of the active column is deposited in the bottom active layer, so the
full incoming shortwave flux is conserved.

When `PenetratingShortwaveTendency.Enable` is true,
`SfcTracerForcingTendencyEnable` must also be enabled. Shortwave is removed
from the surface-only heat sum in `SfcTracerForcingOnCell`; the other surface
heat, freshwater enthalpy, and salt fluxes remain controlled by
`SfcTracerForcingTendencyEnable`.

## Surface tracer restoring

Surface tracer restoring applies a piston-velocity tendency, or damping, at the ocean
surface for selected tracers. This is implemented to mitigate drifts in chosen tracers
(most often salinity) by nudging the model's simulated tracer values towards observed climatological values.
This process prevents oceanic regimes from shifting away from reality due to errors in surface freshwater
forcing (in the case of salinity restoring). Currently, it is applied everywhere when enabled.

### Surface tracer restoring configuration

Surface tracer restoring is controlled by two configuration blocks:

```yaml
Omega:
  SurfaceRestoring:
    TracersToRestore: [Temperature, Salinity]
    PistonVelocity: 1.585e-5

  Tendencies:
    SurfaceTracerRestoringEnable: true
```

- `TracersToRestore`: list of tracer names that restoring is applied to
- `PistonVelocity`: restoring rate coefficient
- `SurfaceTracerRestoringEnable`: switch to enable surface tracer restoring

When restoring is enabled, Omega resolves `TracersToRestore` into an internal
list of tracer IDs and applies restoring only to tracers in that list.

### Restoring target fields

Surface restoring uses auxiliary fields:

- `TracersMonthlySurfClimoCell`: restoring target climatological values
- `SurfTracerRestoringDiffsCell`: computed target-minus-state differences

The restoring tendency is computed at the surface layer only and is limited by
the configured `PistonVelocity` and target-minus-state difference.

## Notes

- If a tracer is not listed in `TracersToRestore`, no restoring tendency is
  applied to that tracer.
- If surface restoring is enabled but no tracer IDs are available at tendency
  compute-time, Omega aborts with an error.
