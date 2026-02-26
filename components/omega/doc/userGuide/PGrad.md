(omega-user-pgrad)=

# Pressure Gradient

The pressure gradient term in the momentum equation represents the force per unit
mass due to horizontal variations in pressure and geopotential. This term is
essential for capturing the dynamics of
ocean circulation, including both barotropic and baroclinic motions.

## Physical Background

In the layered non-Boussinesq momentum equation solved in Omega, the pressure
gradient tendency for each edge and layer includes three contributions:

1. **Montgomery potential gradient**: The horizontal gradient of the Montgomery
   potential ($\alpha p + g z$), averaged across the top and bottom interfaces of
   each layer. The Montgomery potential combines the pressure gradient and the
   geopotential, and its gradient along coordinate surfaces accounts for both the
   direct pressure force and the effect of tilted layer interfaces that arise when
   using a general vertical coordinate.

2. **Specific volume correction**: A correction term proportional to the gradient
   of specific volume (inverse density) at each edge. This term ensures that
   horizontal density variations between the two cells sharing an edge are properly
   represented in the pressure gradient force.

3. **External geopotential forcing**: Contributions from the tidal potential and
   the self-attraction and loading (SAL) terms. These represent gravitational
   forcing from astronomical tides and the deformation of the solid Earth and ocean
   surface in response to the ocean mass distribution. These terms are currently set
   to zero and will be provided by a future tidal forcing module.

## Configuration Options

The pressure gradient method is configured in the input YAML file under the
`PressureGrad` section:

```yaml
PressureGrad:
   PressureGradType: 'centered'
```

### Available Methods

**Centered Difference** (`'centered'` or `'Centered'`)
- Computes the pressure gradient using a centered finite-difference approximation
  of the Montgomery potential gradient and specific volume correction
- Suitable for global ocean simulations without ice shelf cavities
- Default and currently the only fully implemented option

**High-Order** (`'HighOrder1'`)
- Placeholder for a future high-order pressure gradient method based on volume
  integral formulations
- Intended for simulations with ice shelf cavities and steep bathymetry where the
  centered scheme may be inaccurate
- Not yet implemented; selecting this option produces zero pressure gradient tendency

## Dependencies

The pressure gradient calculation requires the following Omega components to be
initialized first:

- [**Horizontal Mesh**](omega-user-horz-mesh): provides mesh geometry including
  distances between cell centers and edge connectivity
- [**Vertical Coordinate**](omega-user-vert-coord): provides pressure at layer
  mid-points and interfaces, interface heights ($z$), and geopotential
- [**Equation of State**](omega-user-eos): provides the specific volume field
- [**Ocean State**](omega-user-state): provides the current layer thicknesses
