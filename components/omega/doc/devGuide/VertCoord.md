(omega-dev-vert-coord)=

## Vertical Coordinate

### Initialization

The default `VertCoord` instance is created by the `init` method, which assumes that `Decomp` has already been initialized.
```
Decomp::init();
VertCoord::init();
```
The `init` method accepts two optional arguments with default values:
```
VertCoord::init(const bool ReadStream = true, const int NVertLayers = 0)
```
These arguments provide flexibility, particularly for unit testing. When `init(false)` is used, the method skips reading
the `InitialVertCoord` stream. This is useful for unit tests that rely on meshes lacking the fields required by that stream.
In this case, the min/max layer arrays are initialized with default values based on the number of vertical layers, and
`bottomDepth` remains uninitialized. Some tests may utilize a specific number of vertical layers that differs from what is
defined in the mesh file. To handle this, the `VertCoord` can be initialized with `init(false, LocNVertLayers)` to explicitly
set the number of vertical layers. An argument for `ReadStream` must be provided in order to explicitly set `NVertLayers`.
For example, `init(42)` is invalid, `init(false, 42)` must be called instead.

The default instance can be retrieved by:
```
auto *DefVertCoord = VertCoord::getDefault();
```

Additional instances can be created by calling the `create` method.
```
VertCoord *VertCoord::create(const std::string &Name,      ///< [in] Name for vertical coordinate
                             const Decomp *MeshDecomp,     ///< [in] Decomp for mesh
                             Config *Options,              ///< [in] Confiuration options
                             const bool ReadStream = true, ///< [in] optional logical to read stream
                             const int NVertLayers = 0     ///< [in] optional int to set vertical dim
)
```
This calls the constructor and places the instance in a map that can be used to retrieve the instance by name:
```
auto *DefVertCoord = VertCoord::get('Name');
```
The constructor stores some information from `Decomp`, allocates the primary member variables that are computed by methods of the `VertCoord` class during a model timestep.
Host mirror copies are also created, with variables appended with `H`.
It also reads in some additional mesh information and computes the information describing the min/max location of the active vertical layers.
Finally, it registers variables with `IOStreams` so they can be requested as output.

### Variables

A list of member variables along with their types and dimension sizes is below:

| Variable Name | Type | Dimensions |
| ------------- | ---- | ---------- |
| PressureInterface | Real | NCellsSize, NVertLayersP1 |
| PressureMid | Real | NCellsSize, NVertLayers |
| ZInterface | Real | NCellsSize, NVertLayersP1|
| ZMid | Real | NCellsSize, NVertLayers |
| GeopotentialMid | Real | NCellsSize, NVertLayers |
| LayerThicknessTarget | Real | NCellsSize, NVertLayers|
| MinLayerCell | Integer | NCellsSize |
| MaxLayerCell | Integer | NCellsSize |
| MinLayerEdgeTop | Integer| NEdgesSize |
| MaxLayerEdgeTop | Integer | NEdgesSize |
| MinLayerEdgeBot | Integer | NEdgesSize |
| MaxLayerEdgeBot | Integer | NEdgesSize |
| MinLayerVertexTop | Integer | NVerticesSize |
| MaxLayerVertexTop | Integer | NVerticesSize |
| MinLayerVertexBot | Integer | NVerticesSize |
| MaxLayerVertexBot | Integer | NVerticesSize |
| VertCoordMovementWeights | Real | NCellsSize, NVertLayers |
| RefPseudoThickness | Real | NCellsSize, NVertLayers |
| BottomDepth | Real | NCellsSize |

### Removal

`VertCoord` instances can be removed by name:
```
VertCoord.erase("Name");
```
or all instances can be destroyed by calling:
```
VertCoord.clear();
```

### Use of hierarchical parallelism

The methods `computePressure` and `computeZHeight` are similar in that they use hierarchical parallelism to split the work for horizontal cells over teams of threads, with a `parallel_for`.
This is done with a `TeamPolicy`:
```c++
const auto Policy = TeamPolicy(NCellsAll, OMEGA_TEAMSIZE, 1);
```
The `parallel_for` is then called with this policy:
```c++
Kokkos::parallel_for("loopName", Policy, KOKKOS_LAMBDA(const TeamMember &Member) {
       const I4 ICell = Member.league_rank();
     ...
}
```
The cumulative sum in the vertical is computed among threads (in parallel) within in the `parallel_for` using  a `parallel_scan`.
The `parallel_scan` is called inside the `parallel_for`:
```c++
Range = KMax - KMin + 1;
Kokkos::parallel_scan(
    TeamThreadRange(Member, Range),
    [=](int K, Real &Accum, bool IsFinal) {
    ...
}
```
where `KMax` and `KMin` are the maximum and minimum active vertical layer indices.
In `computeTargetThickness` an outer loop divides the horizontal cells over teams of threads as above, however, there is a nested `parallel_reduce` that computes the column sum of the reference layer thicknesses times the vertical coordinate movement weights among threads in a team:
```c++
Real SumWh 0= 0;
Kokkos::parallel_reduce(
    Kokkos::TeamThreadRange(Member, KMin, KMax + 1),
    [=](const int K, Real &LocalWh) {
       LocalWh += VertCoordMovementWeights(ICell, K) *
                  RefLayerThickness(ICell, K);
    },
    SumWh);
```
also, the vertical computation of the target thicknesses are computed in a nested `parallel_for`:
```c++
Kokkos::parallel_for(
    Kokkos::TeamThreadRange(Member, NChunks), [=](const int KChunk) {
   ...
}
```
This `parallel_for` iterates over vertical chunks to facilitate vectorization on CPUs within an inner `for` loop over the vector length.
The vector length on GPUs is set to 1 to maximize parallelism.
The `computeGeopotential` method uses hierarchical parallelism in a very similar way to `computeTargetThickness`, except that it doesn't require a column sum.
It has an outer `parallel_for` that splits horizontal cells into teams and an inner `parallel_for` that does vertical computations in chunks.
