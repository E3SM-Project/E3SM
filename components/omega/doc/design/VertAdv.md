(omega-design-vert-adv)=
# Vertical Advection

## 1 Overview

The Vertical Advection module is responsible for the vertical transport of mass,
momentum, and tracers within Omega. This implementation is consistent with
Omega's non-Boussinesq governing equations, where the prognostic vertical
coordinate is a pseudo-height, effectively a normalized pressure coordinate. The
methods compute the tendencies associated with the vertical advection terms in
the governing equations, and are designed to allow for flexibility in numerical
scheme configuration.

## 2 Requirements

### 2.1 Requirement: Compute transport pseudo-velocity

Compute the transport pseudo-velocity through vertical layer interfaces based on
the divergence of the horizontal velocity field. Allow for possible movement of
layer interfaces and tilt of layer interfaces

### 2.2 Requirement: Advection scheme flexibility

Support for multiple possible reconstruction schemes, allowing users to select
among different orders of accuracy and upwinding at runtime through
configuration options. Standard and monotonic flux-corrected transport (FCT)
advection algorithms are available for the tracer tendency calculations.

### 2.3 Requirement: Compute tendencies for mass, velocity, and tracers

Separate methods are needed to compute the vertical-advection tendencies
for mass, velocity, and tracers. Mass (`LayerThickness`) and horizontal
velocity (`NormalVelocity`) are stored as 2-D arrays. `LayerThickness` is
cell-centered, while `NormalVelocity` is defined on edges. These tendencies
use a second-order flux formulation that ensures strict conservation and energy
consistency with the continuity equation. `Tracers` is a 3-D, cell-centered
array. The tracer tendencies are computed using a configurable flux-form scheme
that supports multiple interface reconstruction options.

### 2.4 Requirement: Modular design

The vertical advection module must consist of classes and methods that are easily
incorporated into different time-integration schemes and maintain compatibility
with the Omega driver framework.

### 2.5 Desired: Monotonicity checks

The vertical advection module should include optional monotonicity checks to
ensure that reconstructed tracer fields remain physically bounded and free of
spurious oscillations when using the FCT scheme.

## 3 Algorithmic Formulation

The [Omega V1 Governing Equations](omega-design-governing-eqns-omega1) describe
the evolution of mass, velocity, and tracers, including the contributions
from vertical advection.

### 3.1 Pseudo-velocity

The vertical pseudo-velocity is derived from the divergence of the horizontal
velocity $d\tilde{w} / d\tilde{z} = - \nabla \cdot {\bf u}.$ In
discrete form, the pseudo-velocity at the top interface of layer $k$ in cell
$i$ is computed from the thickness-weighted divergence of horizontal velocity:

$$
\tilde{w}_{i,k-1/2} = \tilde{w}_{i,k+1/2} +
\frac{1}{A_i} \sum_{e\in EC(i)} n_{e,i} [\tilde{h}_{i,k}]_{e} u_{e,k} l_e,
$$

using the TRiSK formulation of the discrete divergence operator (see
[V0 design document](omega-design-shallow-water-omega0)). The boundary
conditions at the top and bottom interfaces set the pseudo-velocity to zero.

The vertical pseudo-velocity $\tilde{w}$ is represented in the code by
`VerticalVelocity`. The equations below use the transport pseudo-velocity
$\tilde{W}_{tr} which incorporates optional corrections on the pseudo-velocity
arising from interface motion, contributions from the horizontal velocity
through tilted interfaces, and surface source terms. This is represented in the
code by `TotalVerticalVelocity`.

### 3.2 Mass (Thickness)

The time evolution of the pseudo-thickness (represented in the code as
`LayerThickness`) due to vertical transport is

$$
\tilde{h}^{n+1}_{i,k} = \tilde{h}^{n}_{i,k} +
\Delta t \left( \left[\tilde{W}_{tr}\right]_{k+1/2} -
\left[\tilde{W}_{tr}\right]_{k-1/2} \right).
$$

This second-order form ensures exact conservation of total column thickness and
numerical stability.

### 3.3 Velocity

The velocity equation governs the evolution of the horizontal velocity field
(`NormalVelocity`),

$$
u^{n+1}_{e,k} = u^{n}_{e,k} +
\frac{\Delta t}{\left[\tilde{h}_{i,k} \right]_{e}} \left( \left[U \tilde{W}_{tr} \right]_{e,k+1/2} -
\left[U \tilde{W}_{tr} \right]_{e,k-1/2} \right),
$$

where $U = u - \tilde{u}$ removes the component of the horizontal velocity normal to
a tilted interface. Square brackets $[\,]$ denote interpolation of quantities
from cell centers to edges or from the middle of a layer to an interface.

### 3.4 Tracers

The contribution of vertical advection to tracer transport is expressed as

$$
\tilde{h}^{n+1}_{i,k} \tilde{\phi}^{n+1}_{i,k} = \tilde{h}^{n}_{i,k} \tilde{\phi}^{n}_{i,k} +
\Delta t \left( \left[ \phi \tilde{W}_{tr} \right]_{k+1/2} -
\left[ \phi \tilde{W}_{tr} \right]_{k-1/2} \right).
$$

A variety of options are available for the scheme used to reconstruct the tracer
values at interfaces $[\phi]_{k\pm1/2}$ (order of accuracy, standard advection vs.
flux-corrected transport) which can be selected at runtime via configuration.

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters

An `enum class` is used for the selection of the reconstruction scheme for the
fluxes at layer interfaces:

```c++
enum class VertFluxOption {
    Second, /// 2nd-order centered
    Third,  /// 3rd-order upwind
    Fourth  /// 4th-order centered
};
```

Another `enum class` specifies the tracer advection scheme:

```c++
enum class VertAdvOption {
    Standard, /// standard non-monotonic scheme
    FCT,      /// flux-corrected transport scheme
};
```

#### 4.1.2 Class/structs/data types

The `VertAdv` class provides the core functionality for vertical advection.
It defines the data structures, configuration parameters, and compute methods
needed for performing advective vertical transport.

```c++
class VertAdv{
   // Logicals to enable/disable advection of specific fields
   bool ThickVertAdvEnabled;
   bool VelVertAdvEnabled;
   bool TracerVertAdvEnabled;

   // Enum options
   VertFluxOption VertFluxChoice;
   VertAdvOption VertAdvChoice;

   // Core data arrays
   Array2DReal VerticalVelocity;      /// pseudo-velocity through top of cell
   Array2DReal TotalVerticalVelocity; /// transport velocity through top of cell
   Array3DReal VerticalFlux;          /// tracer fluxes at vertical interfaces
   Array3DReal LowOrderVerticalFlux;  /// low order fluxes for FCT

   // Compute methods
   computeVerticalVelocity(...);
   computeThicknessVAdvTend(...);
   computeVelocityVadvTend(...);
   computeTracerVadvTend(...);
   computeVerticalFluxes(...);
   computeStdVAdvTend(...);
   computeFCTVAdvTend(...);

   // Mesh dimensions
   I4 NVertLayers;
   I4 NVertLayersP1;
   I4 NCellsOwned;
   I4 NCellsAll;
   I4 NCellsSize;
   I4 NEdgesOwned;
   I4 NEdgesAll;
   I4 NEdgesSize;

   // Vertical loop bounds from VertCoord
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;

   // Arrays from HorzMesh
   Array2DI4 CellsOnEdge;
   Array2DI4 EdgesOnCell;
   Array1DI4 NEdgesOnCell;
   Array1DReal AreaCell;
   Array1DReal DvEdge;
   Array2DReal EdgeSignOnCell;
};
```

### 4.2 Methods

#### 4.2.1 Creation, Destruction, Retrieval

Like other Omega components, instances of the `VertAdv` class will be stored in
a static map and a `DefVAdv` pointer will point to the default instance.

```c++
VertAdv *VertAdv::DefVAdv;
std::map<std::string, std::unique_ptr<VertAdv>> VertAdv::AllVAdv
```

Methods are required for construction, destruction, and retrieval. After
construction, object pointers are inserted in the map; upon destruction, they
are removed. Utility functions provide access to the default `VertAdv` instance
and any non-default instances.

An initialization method ``init`` is responsible for reading configuration
parameters and creating the default `VertAdv` object.

```c++
void VertAdv::init();
```

#### 4.2.2 Pseudo-velocity computation

The quantities required by this method are either stored locally in the
`VertAdv` object, or are fetched from  the `State` and `AuxState` objects.

This method computes the vertical pseudo-velocity based on the divergence of
the horizontal velocity field. A `parallelForOuter` construct iterates over the
cells in the horizontal mesh. Within each cell, a `parallelForInner` iterates
through the active layers to compute the divergence of horizontal velocity,
storing the result in scratch space. After the divergence is computed, a
`parallelScanInner` performs a prefix sum to obtain the vertical pseudo-velocity
at each layer interface. These pseudo-velocities are stored in the
`VerticalVelocity` member array. The transport pseudo-velocity
`TotalVerticalVelocity` is then computed from the `VerticalVelocity` and
any desired corrections. `TotalVerticalVelocity` is used by the tendency
methods below.


```c++
void VertAdv::computeVerticalVelocity(const Array2DReal &NormalVelocity,
    const Array2DReal &FluxLayerThickEdge, TimeInterval Dt)
```

#### 4.2.3 Tendency computations

The `computeThicknessVAdvTend`, `computeVelocityVadvTend`, and
`computeTracerVadvTend` methods are designed to be called by the
`Tendencies` class methods.

The `computeThicknessVAdvTend` method computes the vertical advection tendency
for pseudo-thickness (`LayerThickness`). The parallel execution of the
computation is straightforward.

```c++
void VertAdv::computeThicknessVAdvTend(const Array2DReal &Tend) {

   if (!ThickVertAdvEnabled) return;

   OmegaScope(LocTotVertVel, TotalVerticalVelocity);

   parallelForOuter(...,
      parallelForInner(...,
         Tend(ICell, K) += LocTotVertVel(ICell, K + 1) -
                           LocTotVertVel(ICell, K);
      );
   );
}
```

The `computeVelocityVadvTend` method computes the vertical advection tendency
for horizontal velocity (`NormalVelocity`). The method requires the
`NormalVelocity` field from the `State` object and the `FluxLayerThickEdge`
field from the `AuxState` object as input. A `parallelForOuter` loop iterates
over the edges in the horizontal mesh, and a `parallelForInner` loop iterates
over the active layer interfaces. At each interface, `NormalVelocity` and
`FluxLayerThickEdge` are used to compute $du/dz$, and `TotalVerticalVelocity`
is interpolated from cell centers to edges. The product of these, $W du/dz$ is
averaged from interfaces to the middle of each layer and accumulated into the
velocity tendency.

```c++
void VertAdv::computeVelocityVAdvTend(const Array2DReal &Tend,
    const Array2DReal &NormalVelocity, const Array2DReal &FluxLayerThickEdge)
```

The `computeTracerVAdvTend` method serves as an interface that calls
`computeVerticalFluxes` to compute the needed tracer fluxes based on specified
configuration options and then dispatches to either the standard advection
algorithm, or the FCT scheme depending on runtime configuration.

```c++
void VertAdv::computeTracerVadvTend(const Array3DReal &Tend,
    const Array3DReal &Tracers)
```

The `computeStdVAdvTend` method performs the standard advection update by
looping over the cells and interfaces using nested `parallelForOuter` and
`parallelForInner` constructs.

The `computeFCTVAdvTend` method implements the flux-corrected transport scheme
developed by [Zalesak 1979](https://www.sciencedirect.com/science/article/pii/0021999179900512).
We use by default the FCT formulation of
[Skamarack & Gassmann 2011](https://journals.ametsoc.org/view/journals/mwre/139/9/mwr-d-10-05056.1.xml).
At each interface, both a high-order flux (chosen via configuration) and a
low-order flux (1st-order upwind) are computed in `computeVerticalFluxes`.
These fluxes are then blended and limited to ensure monotone transport.

## 5 Verification and Testing

Unit tests are required for each compute method to verify correct operation and
numerical consistency. In addition, convergence tests with the
[merry-go-round](https://docs.e3sm.org/polaris/main/users_guide/ocean/tasks/merry_go_round.html)
test case should confirm the expected order of accuracy of the tracer fluxes.
