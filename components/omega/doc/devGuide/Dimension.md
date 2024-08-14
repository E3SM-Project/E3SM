(omega-dev-dimension)=

## Dimension

When performing IO for fields in Omega, some metadata associated with each
dimension of multi-dimensional arrays are required. The Dimension class is
a container that stores this dimension information and a Dimension instance
must be created for each dimension in Omega. Most of these will be created
by the Mesh or Decomposition classes, but if developers introduce a new
dimension, the dimension must be created.

For dimensions that are distributed across the parallel machine, the creation
must provide length and an offset array that describes the global index offset
for each local point along the dimension. Typically this is the zero-based
ID of the location in the global mesh. A negative value is used for any halo
regions that will not participate in IO. For example, the NCells dimension
is created using:
```c++
   Decomp *DefDecomp = Decomp::getDefault();

   I4 NCellsSize   = DefDecomp->NCellsSize; // Local dim length
   I4 NCellsGlobal = DefDecomp->NCellsGlobal; // Global dim length
   I4 NCellsOwned  = DefDecomp->NCellsOwned;

   HostArray1DI4 CellOffset("NCellsOffset", NCellsSize);

   for (int N = 0; N < NCellsSize; ++N) {
      if (N < NCellsOwned) {
         CellOffset(N) = DefDecomp->CellIDH(N) - 1; // Offset is zero-based
      } else {
         CellOffset(N) = -1; // Denotes halo cells not to be used in IO
      }
   }

   std::shared_ptr<Dimension> CellDim =
      Dimension::create("NCells", NCellsGlobal, NCellsSize, CellOffset);
```

For dimensions that are local (eg the vertical dimension or MaxEdgeOnCell),
only the name and length are needed for creation using:
```c++
   std::shared_ptr<Dimension> LocalDim =
      Dimension::create("MyDimName", LengthOfDim);
```
The offset array is computed internally for these dimensions and is simply
the array index along the dimension.

Each Dimension stores the global and local length of the dimension as well
as the offset array. It also stores a flag to denote whether the dimension
is distributed or not. This information is used primarily by the IOStreams
class to write dimension metadata and define offsets for parallel IO. It is
not expected to be used elsewhere in Omega.

Dimension information is available through retrieval functions. These can be
retrieved by dimension name as in:
```c++
   I4 LengthGlob  = Dimension::getDimLengthGlobal("MyDimName");
   I4 LengthLoc   = Dimension::getDimLengthLocal("MyDimName");
   bool DistrbDim = Dimension::isDistributedDim("MyDimName");
   HostArray1DI4 MyOffset = Dimension::getDimOffset("MyDimName");
```
or it can be retrieved from an instance directly using member retrieval
functions:
```c++
   std::shared_ptr<Dimension> MyDim = Dimension::get("MyDimName");
   I4 LengthGlob = MyDim->getLengthGlobal();
   I4 LengthLoc  = MyDim->getLengthLocal();
   bool DistrbDim = MyDim->isDistributed();
   HostArray1DI4 MyOffset = MyDim->getOffset();
```
An iterator is also provided to enable looping over all defined dimensions:
```c++
   for (auto Iter = Dimension::begin(); Iter != Dimension::end(); ++Iter) {
      std::string MyDimName = Iter->first;
      std::shared_ptr<Dimension> ThisDim = Iter->second;
      I4 LengthLoc = ThisDim->getLengthLocal();
      // other retrievals, etc
   }
```

As in other classes, a Dimension can be removed using:
```c++
   Dimension::destroy("MyDimName");
```
and all dimensions should be removed during Omega finalization before exiting
the Kokkos environment using:
```c++
   Dimension::clear();
```

Examples of the use of Dimensions can be found in the Dimension unit test.
