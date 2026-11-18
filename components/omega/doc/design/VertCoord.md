(omega-design-vert-coord)=
# Vertical Coordinate

## 1 Overview
The vertical coordinate module will be responsible for computing and storing information related to the vertical mesh in Omega.
Since Omega will be a non-Boussinesq model, the vertical independent variable will be pressure-based (expressed in terms of the pseudo-height) as opposed to $z$-based.
Similar to MPAS-Ocean's support for tilted height ($z^\star$) coordinates, Omega V1 will support a general vertical coordinate, where pseudo thickness $\tilde{h}$ varies in proportion to the total pseudo thickness of the water column.
The prognostic variable in the layered continuity equation is pseudo thickness, which can be used to calculate the total pressure at a given vertical layer interface.
The geometric height of a given layer interface is still necessary in the model to compute sea level, the geopotential, and derivatives with respect to $z$, all of which affect the dynamics.
The geometric thickness of a layer ($\Delta z$) can be found using the pseudo thickness, density, reference density, and known bottom elevation (positive up) relative to the reference geoid.
The vertical coordinate module will also serve as a container for information related to the variable number of active vertical layers for a given ocean column (due to variations in bottom elevation and ice shelf cavities).


## 2 Requirements

### 2.1 Requirement: Compute and store the total pressure at a given layer interface

The total pressure at each layer interface will be computed and stored within the vertical coordinate module.
This involves a surface-down summation of the pseudo thickness variable times the gravitational acceleration, $g$, and the reference density, $\rho_0$, starting with the surface pressure.
The pressure variable will be needed in the computation of the TEOS-10 equation of state and the pressure gradient.
The variable will be registered with the `IOStreams` framework to enable it to be output during runtime.

### 2.2 Requirement: Compute and store the $z$ location of a given layer interface

The $z$ location at each layer interface will be computed and stored within the vertical coordinate module.
In the model, $z$ is defined as being positive up from the geoid.
The $z$ location is calculated in a bottom-up summation of the pseudo thickness times the specific volume and reference density, starting with the known bottom elevation.
The $z$ location of the interfaces is needed to compute the geopotential gradient and the sea level.
The variable will be registered with the `IOStreams` framework to enable it to be output during runtime.

### 2.3 Requirement: Compute and store the vertical layer bounds for each horizontal cell

All tendency terms require a variable vertical loop range to account for variation in bottom elevation and ice shelf cavities.
These bounds will be computed on construction in the vertical coordinate module and used for setting the loop bounds for active layers in calculations over the vertical.

### 2.4 Requirement: Compute and store the geopotential

The geopotential will be computed and stored within the vertical coordinate module.
Initially this will be calculated by summing together the $z$ location times $g$.
In the future, tidal potential and self attraction and loading will also have optional (default off) contributions to the geopotential.


### 2.5 Requirement: Compute and store the desired thickness used in the $p^\star$ vertical coordinate

The $p^\star$ vertical coordinate specifies a desired layer thickness that stretches the initial thicknesses by the bottom pressure anomaly.
This involves adding a vertically-weighted bottom pressure perturbation to the reference (initial) thickness.
This desired thickness will be computed and stored within the vertical coordinate module.
The vertical advection algorithm will use this desired thickness to compute the vertical transport.

### 2.6 Desired: General Arbitrary Lagrangian Eulerian (ALE) coordinate support

In future versions of Omega, the $p^\star$ coordinate will be extended to a more general ALE coordinate

### 2.7 Desired: Vertical Lagrangian remapping (VLR) support

The vertical coordinate module should retain flexibility to support VLR in future versions of Omega.

## 3 Algorithmic Formulation
In the layered non-Boussinesq equations solved in Omega (see [V0 governing equation document](omega-design-governing-eqns-omega1) for details), the prognostic variable for cell $i$ and layer $k$ is the pseudo thickness, $\tilde{h}_{i,k}$, so that the geometric thickness (in meters) is a diagnostic variable defined as:

$$ \Delta z_{i,k} = \rho_0 \alpha_{i,k} \tilde{h}_{i,k}. $$

The pressure at vertical cell interfaces is the found by summing the pseudo thicknesses:

$$  p_{i,k+1/2} = p_{i}^{surf} + g\rho_0 \sum_{k^\prime=1}^k \tilde{h}_{i,k^\prime}. $$

and at the midpoint by

$$ p_{i,k} = p_{i}^{surf} + g\rho_0 \sum_{k'=1}^{k-1} \tilde{h}_{i,k'} + \frac{1}{2} g\rho_0 \tilde{h}_{i,k} $$

The $z$ location of cell interfaces is found by summing the pseudo thicknesses from the bottom multiplied by the specific volume:

$$ z_{i,k+1/2} = z_i^{floor} + \rho_0 \sum_{k^\prime=k}^{K_{max}} \alpha_{i,k^\prime} \tilde{h}_{i,k^\prime}, $$

where $z_i^{floor}$ is the (positive-up) bottom elevation.
The $z$ location of a layer midpoint is given by:

$$ z_{i,k} = z_i^{floor} + \frac{1}{2} \rho_0\alpha_{i,k} \tilde{h}_{i,k} + \rho_0\sum_{k^\prime= k+1}^{K_{max}} \alpha_{i,k^\prime} \tilde{h}_{i,k^\prime}. $$

Initially, the geopotential is the sum of the $z$ height times $g$.
In the future, it will include contributions from the tidal potential ($\Phi_{tp}$) and self attraction and loading ($\Phi_{SAL}$):

$$ \Phi_{i,k} = \left( gz_{i,k} + \Phi_{tp} + \Phi_{SAL} \right). $$

The desired pseudo thickness used for the $p^\star$ coordinate is:

$$\tilde{h}_k^{p^\star} = \tilde{h}_k^{ref} + \left(\frac{p_B-p_{surf}}{g\rho_0} - \sum_{k^\prime=1}^K \tilde{h}_{k^\prime}^{ref} \right)\frac{W_k\tilde{h}_k^{ref}}{\sum_{k^\prime=1}^K W_{k^\prime}\tilde{h}_{k^\prime}^{ref}}, $$

where $\tilde{h}_k^{ref}$ is the initial reference pseudo thickness and the weights $W_k$ determine how pressure perturbations are distributed amongst the layers.
Setting all $W_k$ values to a constant, corresponds to uniform stretching of the layers.
Different weight values can be chosen such that perturbations are distributed unequally.


## 4 Design

### 4.1 Data types and parameters

The `VertCoord` class will be used to compute pressure, $z$ height, and geopotential at each layer.
It will also contain the vertical loop bounds for cells, edges, and vertices.
```c++
class VertCoord {
    public:
        I4 NVertLevels;
        I4 NVertLevelsP1;

        // Variables computed
        Array2DReal PressureInterface;
        Array2DReal PressureMid;
        Array2DReal ZInterface;
        Array2DReal ZMid;
        Array2DReal GeopotentialMid;
        Array2DReal LayerThicknessTarget;

        // Vertical loop bounds (computed on construction)
        Array1DI4 MinLevelCell;
        Array1DI4 MaxLevelCell;
        Array1DI4 MinLevelEdgeTop;
        Array1DI4 MaxLevelEdgeTop;
        Array1DI4 MinLevelEdgeBot;
        Array1DI4 MaxLevelEdgeBot;
        Array1DI4 MinLevelVertexTop;
        Array1DI4 MaxLevelVertexTop;
        Array1DI4 MinLevelVertexBot;
        Array1DI4 MaxLevelVertexBot;

        // Variables read in from vert coord stream
        /// p star coordinate variables
        Array2DReal VertCoordMovementWeights;
        Array2DReal RefLayerThickness;

        /// Variables read in from mesh file
        Array2DReal BottomDepth;
    private:

        // Variables from HorzMesh
        I4 NCells;
        I4 NEdges;
        I4 NVertices;
        Array2DReal CellsOnEdge;
        Array2DReal CellsOnVertex;

        static VertCoord *DefaultVertCoord;
        static std::map<std::string, std::unique_ptr<VertCoord>> AllVertCoords;
}
```

### 4.2 Methods

```c++
class VertCoord{
    public:
        static VertCoord *create();
        void computePressure();
        void computeZHeight();
        void computeGeopotential();
        void computeTargetThickness();
    private:
        void minMaxLevelEdge();
        void minMaxLevelVertex();

}
```
#### 4.2.1 Creation

The constructor will be responsible for storing any static mesh information as private variables and handling the selection any user-specified options.
It will also compute the vertical loop bounds on edges and vertices.
```c++
VertCoord::VertCoord(const HorzMesh *Mesh,
                     int NVertLevels,
                     Config *Options);
```

The create method will take the same arguments as the constructor plus a name.
It calls the constructor to create a new vertical coordinate instance, and put it in the static map of all vertical coordinate objects.
It will return a pointer to the newly created object.
```c++
VertCoord *VertCoord::create(const std::string &Name,
                             const HorzMesh *Mesh,
                             int NVertLevels,
                             Config *Options);
```
A vertical coordinate `IOStream` will be defined to read in the `BottomDepth`, `VertCoordMovementWeights` and `RefLayerThickness` variables.

#### 4.2.2 Initialization

The init method will create the default vertical coordinate and return an error code:
```c++
int VertCoord::init();
```

#### 4.2.3 Retrieval

There will be methods for getting the default and non-default vertical coordinate instances:
```c++
VertCoord *VertCoord::getDefault();
VertCoord *VertCoord::get(const std::string &Name);
```

#### 4.2.4 Computation

The public `computePressure` will sum the pseudo thicknesses times $g$ from the top layer down, starting with the surface pressure.
This will be done with a `parallel_scan` inside a `parallel_for` over the mesh cells using hierarchical parallelism.
```c++
void VertCoord::computePressure(const Array2DReal &PressureInterface,
                                const Array2DReal &PressureMid,
                                const Array2DReal &LayerThickness,
                                const Array2DReal &SurfacePressure) {

}
```

The public `computeZHeight` will sum the pseudo thicknesses times $\alpha$ from the bottom later up, starting with the bottom elevation.
This will be done with a `parallel_scan` inside a `parallel_for` over the mesh cells using hierarchical parallelism.
```c++
void VertCoord::computeZHeight(const Array2DReal &ZInterface,
                               const Array2DReal &ZMid,
                               const Array2DReal &LayerThickness,
                               const Array2DReal &SpecVol,
                               const Array2DReal &BottomDepth) {

}
```

The public `computeGeopotential` will sum together the $z$ height times $g$, the tidal potential, and self attraction and loading:
```c++
void VertCoord::computeGeopotential(const Array2DReal &GeopotentialMid,
                                    const Array2DReal &ZMid,
                                    const Array2DReal &TidalPotential
                                    const Array2DReal &SelfAttractionLoading) {

}
```
Tidal potential forcing and self attraction and loading will be default-off features.
The will be added to (or excluded from) the geopotential based on config flags.

The public `computeTargetThickness` will determine the desired pseudo thickness used for the $p^\star$ vertical coordinate:
```c++
void VertCoord::computeTargetThickness(const Array2DReal &LayerThicknessPStar,
                                       const Array2DReal &VertCoordMovementWeights,
                                       const Array2DReal &RefLayerThickness) {

}
```

The private `minMaxLevel` will determine the various vertical loop bounds on edges and vertices.
It will compute the class member variables: `MinLevelEdgeTop`, `MaxLevelEdgeTop`, `MinLevelEdgeBot`, `MaxLevelEdgeBot`, `MinLevelVertexTop`, `MaxLevelVertexTop`, `MinLevelVertexBot`, and `MaxLevelVertexBot` from `MinLevelCell`, `MaxLevelCell`, `CellsOnEdge`, and `CellsOnVertex`.
```c++
void VertCoord::minMaxLevel( ) {

}
```

#### 4.2.5 Destruction and removal

No operations are needed in the destructor.
The erase method will remove a named vertical coordinate instance, whereas the clear method will remove all of
them.
Both will call the destructor in the process.
```c++
void VertCoord::erase(const std::string &Name);
void VertCoord::clear();
```

## Verification and Testing

### Test:
Unit tests will be used to test each of the computations (computePressure, computeZHeight, computeGeopotential, computePStarThickness) for a given pseudo thickness array. Comparison will be made to a ''truth'' vertical mesh that includes spatially varying `minLevelCell` and `maxLevelCell`.
