(omega-user-frazil)=

# Frazil

This page describes user-facing configuration for frazil tendencies in Omega.
Frazil physics is used to represent the formation and melt of frazil ice within
the ocean water column. It impacts the local layer pseudo-thickness (i.e. mass),
temperature, and salinity tendencies. The vertical sum of the frazil energy, mass of water and mass of salt are passed to the coupler (if coupled) or discarded (in ocean standalone mode).

## Scientific basis
Frazil ice forms in the water column wherever the water becomes colder than its local freezing temperature (aka super-cooled water). The implementation in Omega checks each column for the presence of such super-cooled water. Where it is found, the excess heat deficit is converted into frazil ice: the layer gains heat (and so warms back towards the freezing point), loses mass to the ice, and becomes saltier because frazil has a lower average salinity than seawater. Because frazil ice is buoyant, it is assumed to rise, so an ocean column is processed from the bottom upwards and ice formed at depth is added to a reservoir to interact with the layers above. If one of those layers is warmer than its freezing point, some or all of the rising ice melts there, cooling and freshening that layer. Only the net frazil mass, salt, and energy (i.e. the frazil that survives to the top) are passed to the coupler (or discarded in standalone mode). The optional `DepthLimit` parameter confines this process to the upper part of the column, to avoid instability or unrealistically deep frazil formation. Although the frazil code is designed to bring super-cooled water to its freezing temperature instantly , its implementation inside multi-phase time steppers means super-cooled water is eliminated gradually over several time steps. The details of the salt, mass and energy changes with the frazil formation (and its melt) is dependent on the option chosen (`FixedPropery` vs `Teos`).

## Configuration overview

Frazil behavior is controlled by one enable switch in `Tendencies` and one
`Frazil` configuration block:

```yaml
Omega:
  Tendencies:
    FrazilTendencyEnable: true

  Frazil:
    FrazilType: FixedProperty
    LayerMassFracMax: 0.1
    Phi: 0.75
    DepthLimit: -1.0
    ConservationCheck: false
```

- `Tendencies.FrazilTendencyEnable`
  - enables/disables application of frazil tendency contributions
- `Frazil.FrazilType`
  - selects frazil option
  - supported options in current code: `FixedProperty` and `Teos10`
- `Frazil.LayerMassFracMax`
  - limits per-layer frazil mass/thickness tendency magnitude (applied to formation and melt)
- `Frazil.Phi`
  - liquid-mass fraction of frazil (used by teos frazil formation, or by fixed-property frazil when using porosity rather than constant salinity).
- `Frazil.DepthLimit`
  - limits the depth range over which frazil is computed
  - negative values mean no depth limit
- `Frazil.ConservationCheck`
  - enables a column-level diagnostic conservation check with logging

## Available frazil options

Omega currently includes two active frazil pathways:

- `FixedProperty`
  - freezing is based on the formation of fresh solid ice, to which salt is added (similar the mpas-ocean implementation).
- `Teos10`
  - teos-10-based option using Gibbs SeaWater thermodynamic routines.

Both pathways contribute to:

- pseudo-thickness tendency
- temperature tracer tendency
- salinity tracer tendency

## Notes

- Frazil tendencies are applied through the `Tendencies` tracer-step workflow.
- The frazil tendency hook assumes tracer names include `Temperature` and
  `Salinity`.
- For implementation and algorithm details, see
  [Developer Frazil Guide](../devGuide/Frazil.md).
