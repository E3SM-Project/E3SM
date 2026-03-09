//===-- base/HorzMesh.cpp - horizontal mesh methods -------------*- C++ -*-===//
//
// The mesh (Mesh) class initializes a local mesh domain based on a given
// domain decomposition. It retrieves the mesh count and connectivity
// information from the Decomp object and reads in all other mesh variables
// from the mesh file. It also manages the device copies of the mesh data.
// It is meant to provide a container for passing mesh variables throughout
// the OMEGA tendency computation routines.
//
//===----------------------------------------------------------------------===//

#include "HorzMesh.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "IO.h"
#include "Logging.h"
#include "OmegaKokkos.h"

namespace OMEGA {

// create the static class members
HorzMesh *HorzMesh::DefaultHorzMesh = nullptr;
std::map<std::string, std::unique_ptr<HorzMesh>> HorzMesh::AllHorzMeshes;

//------------------------------------------------------------------------------
// Initialize the mesh. Assumes that Decomp and VertCoord have already been
// initialized.

void HorzMesh::init() {

   Error Err; // default successful error code

   // Retrieve the default decomposition
   Decomp *DefDecomp = Decomp::getDefault();

   // Create the default mesh and set pointer to it
   HorzMesh::DefaultHorzMesh = create("Default", DefDecomp);
}

//------------------------------------------------------------------------------
// Construct a new local mesh given a decomposition

HorzMesh::HorzMesh(const std::string &Name, //< [in] Name for new mesh
                   Decomp *MeshDecomp       //< [in] Decomp for the new mesh
                   )
    : CellID(MeshDecomp->CellID) ///< global cell ID for each local cell
{
   MeshName = Name;

   // Retrieve mesh files name from Decomp
   MeshFileName = MeshDecomp->MeshFileName;

   // Retrieve mesh cell/edge/vertex totals from Decomp
   NCellsHalo   = MeshDecomp->NCellsHalo;
   NCellsHaloH  = MeshDecomp->NCellsHaloH;
   NCellsOwned  = MeshDecomp->NCellsOwned;
   NCellsAll    = MeshDecomp->NCellsAll;
   NCellsSize   = MeshDecomp->NCellsSize;
   NCellsGlobal = MeshDecomp->NCellsGlobal;

   NEdgesHalo     = MeshDecomp->NEdgesHalo;
   NEdgesHaloH    = MeshDecomp->NEdgesHaloH;
   NEdgesOwned    = MeshDecomp->NEdgesOwned;
   NEdgesAll      = MeshDecomp->NEdgesAll;
   NEdgesSize     = MeshDecomp->NEdgesSize;
   MaxCellsOnEdge = MeshDecomp->MaxCellsOnEdge;
   MaxEdges       = MeshDecomp->MaxEdges;
   NEdgesGlobal   = MeshDecomp->NEdgesGlobal;

   NVerticesHalo  = MeshDecomp->NVerticesHalo;
   NVerticesHaloH = MeshDecomp->NVerticesHaloH;
   NVerticesOwned = MeshDecomp->NVerticesOwned;
   NVerticesAll   = MeshDecomp->NVerticesAll;
   NVerticesSize  = MeshDecomp->NVerticesSize;
   VertexDegree   = MeshDecomp->VertexDegree;

   // Retrieve connectivity arrays from Decomp
   CellsOnCellH    = MeshDecomp->CellsOnCellH;
   EdgesOnCellH    = MeshDecomp->EdgesOnCellH;
   NEdgesOnCellH   = MeshDecomp->NEdgesOnCellH;
   VerticesOnCellH = MeshDecomp->VerticesOnCellH;
   CellsOnEdgeH    = MeshDecomp->CellsOnEdgeH;
   EdgesOnEdgeH    = MeshDecomp->EdgesOnEdgeH;
   NEdgesOnEdgeH   = MeshDecomp->NEdgesOnEdgeH;
   VerticesOnEdgeH = MeshDecomp->VerticesOnEdgeH;
   CellsOnVertexH  = MeshDecomp->CellsOnVertexH;
   EdgesOnVertexH  = MeshDecomp->EdgesOnVertexH;

   CellsOnCell    = MeshDecomp->CellsOnCell;
   EdgesOnCell    = MeshDecomp->EdgesOnCell;
   NEdgesOnCell   = MeshDecomp->NEdgesOnCell;
   VerticesOnCell = MeshDecomp->VerticesOnCell;
   CellsOnEdge    = MeshDecomp->CellsOnEdge;
   EdgesOnEdge    = MeshDecomp->EdgesOnEdge;
   NEdgesOnEdge   = MeshDecomp->NEdgesOnEdge;
   VerticesOnEdge = MeshDecomp->VerticesOnEdge;
   CellsOnVertex  = MeshDecomp->CellsOnVertex;
   EdgesOnVertex  = MeshDecomp->EdgesOnVertex;

   // Open the mesh file for reading (assume IO has already been initialized)
   IO::openFile(MeshFileID, MeshFileName, IO::ModeRead);

   // Create Omega Dimensions associated with this mesh
   createDimensions(MeshDecomp);

   // Temporary - remove when IOStreams implemented
   // Create the parallel IO decompositions required to read in mesh variables
   initParallelIO(MeshDecomp);

   // Read x/y/z and lon/lat coordinates for cells, edges, and vertices
   readCoordinates();

   // Read the cell-centered bottom depth
   readBottomDepth();

   // Read the mesh areas, lengths, and angles
   readMeasurements();

   // Read the edge mesh weights
   readWeights();

   // Read the Coriolis parameter at the cells, edges, and vertices
   readCoriolis();

   // Destroy the parallel IO decompositions
   finalizeParallelIO();

   // Copy host data to device
   copyToDevice();

   // TODO: add ability to compute (rather than read in)
   // dependent mesh quantities

   // Compute EdgeSignOnCells and EdgeSignOnVertex
   computeEdgeSign();

   // set mesh scaling coefficients
   setMeshScaling();

} // end horizontal mesh constructor

/// Creates a new mesh by calling the constructor and puts it in the
/// AllHorzMeshes map
HorzMesh *HorzMesh::create(const std::string &Name, //< [in] Name for new mesh
                           Decomp *MeshDecomp //< [in] Decomp for the new mesh
) {
   // Check to see if a mesh of the same name already exists and
   // if so, exit with an error
   if (AllHorzMeshes.find(Name) != AllHorzMeshes.end()) {
      LOG_ERROR("Attempted to create a HorzMesh with name {} but a HorzMesh of "
                "that name already exists",
                Name);
      return nullptr;
   }

   // create a new mesh on the heap and put it in a map of
   // unique_ptrs, which will manage its lifetime
   auto *NewHorzMesh = new HorzMesh(Name, MeshDecomp);
   AllHorzMeshes.emplace(Name, NewHorzMesh);

   return NewHorzMesh;
}

//------------------------------------------------------------------------------
// Destroys a local mesh and deallocates all arrays
HorzMesh::~HorzMesh() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end deconstructor

//------------------------------------------------------------------------------
// Removes a mesh from list by name
void HorzMesh::erase(std::string InName // [in] name of mesh to remove
) {

   AllHorzMeshes.erase(InName); // remove the mesh from the list and in
                                // the process, calls the destructor

} // end mesh erase
//------------------------------------------------------------------------------
// Removes all meshes to clean up before exit
void HorzMesh::clear() {

   AllHorzMeshes.clear(); // removes all meshes from the list and in
                          // the porcess, calls the destructors for each

} // end clear

//------------------------------------------------------------------------------
// Create Dimension instances for all dimensions associated with the mesh
void HorzMesh::createDimensions(Decomp *MeshDecomp) {

   // Create non-distributed dimensions
   // If these have already been created (eg not the default mesh), skip
   if (!Dimension::exists("MaxCellsOnEdge"))
      auto MaxCellsOnEdgeDim =
          Dimension::create("MaxCellsOnEdge", MaxCellsOnEdge);
   if (!Dimension::exists("MaxEdges"))
      auto MaxEdgesDim = Dimension::create("MaxEdges", MaxEdges);
   if (!Dimension::exists("VertexDegree"))
      auto VertexDegreeDim = Dimension::create("VertexDegree", VertexDegree);

   // For distributed dimensions we need to compute offsets along each
   // dimension (the global offset at each local point with -1 for non-owned
   // points)

   // Distinguish these dimensions by name if this is not the default mesh
   std::string NameSuffix;
   if (MeshName == "Default") {
      NameSuffix = "";
   } else {
      NameSuffix = MeshName;
   }

   // Create the offset and Dimension for NCells
   std::string DimName = "NCells" + NameSuffix;
   HostArray1DI4 CellOffset("CellOffSet", NCellsSize);
   I4 NCellsGlobal = MeshDecomp->NCellsGlobal;
   // NCellsOwned, NCellsSize already set
   for (int Cell = 0; Cell < NCellsSize; ++Cell) {
      if (Cell < NCellsOwned) {
         CellOffset(Cell) = MeshDecomp->CellIDH(Cell) - 1;
      } else {
         CellOffset(Cell) = -1;
      }
   }
   // Create the distributed dimension
   auto NCellsDim =
       Dimension::create(DimName, NCellsGlobal, NCellsSize, CellOffset);

   // Create the offset and Dimension for NEdges
   DimName = "NEdges" + NameSuffix;
   HostArray1DI4 EdgeOffset("EdgeOffSet", NEdgesSize);
   I4 NEdgesGlobal = MeshDecomp->NEdgesGlobal;
   // NEdgesOwned, NEdgesSize already set
   for (int Edge = 0; Edge < NEdgesSize; ++Edge) {
      if (Edge < NEdgesOwned) {
         EdgeOffset(Edge) = MeshDecomp->EdgeIDH(Edge) - 1;
      } else {
         EdgeOffset(Edge) = -1;
      }
   }
   // Create the distributed dimension
   auto NEdgesDim =
       Dimension::create(DimName, NEdgesGlobal, NEdgesSize, EdgeOffset);

   // Create the offset and Dimension for NVertices
   DimName = "NVertices" + NameSuffix;
   HostArray1DI4 VrtxOffset("VertexOffSet", NVerticesSize);
   I4 NVerticesGlobal = MeshDecomp->NVerticesGlobal;
   // NVerticesOwned, NVerticesSize already set
   for (int Vrtx = 0; Vrtx < NVerticesSize; ++Vrtx) {
      if (Vrtx < NVerticesOwned) {
         VrtxOffset(Vrtx) = MeshDecomp->VertexIDH(Vrtx) - 1;
      } else {
         VrtxOffset(Vrtx) = -1;
      }
   }
   // Create the distributed dimension
   auto NVerticesDim =
       Dimension::create(DimName, NVerticesGlobal, NVerticesSize, VrtxOffset);

} // end createDimensions

//------------------------------------------------------------------------------
// Initialize the parallel IO decompositions for the mesh variables
void HorzMesh::initParallelIO(Decomp *MeshDecomp) {

   I4 NDims             = 1;
   IO::Rearranger Rearr = IO::RearrBox;

   // Create the IO decomp for arrays with (NCells) dimensions
   std::vector<I4> CellDims{MeshDecomp->NCellsGlobal};
   std::vector<I4> CellID(NCellsAll);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      CellID[Cell] = MeshDecomp->CellIDH(Cell) - 1;
   }

   CellDecompR8 = IO::createDecomp(IO::IOTypeR8, NDims, CellDims, NCellsAll,
                                   CellID, Rearr);

   // Create the IO decomp for arrays with (NEdges) dimensions
   std::vector<I4> EdgeDims{MeshDecomp->NEdgesGlobal};
   std::vector<I4> EdgeID(NEdgesAll);
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      EdgeID[Edge] = MeshDecomp->EdgeIDH(Edge) - 1;
   }

   EdgeDecompR8 = IO::createDecomp(IO::IOTypeR8, NDims, EdgeDims, NEdgesAll,
                                   EdgeID, Rearr);

   // Create the IO decomp for arrays with (NVertices) dimensions
   std::vector<I4> VertexDims{MeshDecomp->NVerticesGlobal};
   std::vector<I4> VertexID(NVerticesAll);
   for (int Vertex = 0; Vertex < NVerticesAll; ++Vertex) {
      VertexID[Vertex] = MeshDecomp->VertexIDH(Vertex) - 1;
   }

   VertexDecompR8 = IO::createDecomp(IO::IOTypeR8, NDims, VertexDims,
                                     NVerticesAll, VertexID, Rearr);

   // Create the IO decomp for arrays with (NEdges, 2*MaxEdges) dimensions
   NDims     = 2;
   MaxEdges2 = 2 * MaxEdges;
   std::vector<I4> OnEdgeDims2{MeshDecomp->NEdgesGlobal, MaxEdges2};
   I4 OnEdgeSize2 = NEdgesAll * MaxEdges2;
   std::vector<I4> OnEdgeOffset2(OnEdgeSize2, -1);
   for (int Edge = 0; Edge < NEdgesAll; Edge++) {
      for (int i = 0; i < MaxEdges2; i++) {
         I4 GlobalID = EdgeID[Edge] * MaxEdges2 + i;

         OnEdgeOffset2[Edge * MaxEdges2 + i] = GlobalID;
      }
   }

   OnEdgeDecompR8 = IO::createDecomp(IO::IOTypeR8, NDims, OnEdgeDims2,
                                     OnEdgeSize2, OnEdgeOffset2, Rearr);

   // Create the IO decomp for arrays with (NVertices, VertexDegree) dimensions
   std::vector<I4> OnVertexDims{MeshDecomp->NVerticesGlobal, VertexDegree};
   I4 OnVertexSize = NVerticesAll * VertexDegree;
   std::vector<I4> OnVertexOffset(OnVertexSize, -1);
   for (int Vertex = 0; Vertex < NVerticesAll; Vertex++) {
      for (int i = 0; i < VertexDegree; i++) {
         I4 GlobalID = VertexID[Vertex] * VertexDegree + i;
         OnVertexOffset[Vertex * VertexDegree + i] = GlobalID;
      }
   }

   OnVertexDecompR8 = IO::createDecomp(IO::IOTypeR8, NDims, OnVertexDims,
                                       OnVertexSize, OnVertexOffset, Rearr);

} // end initParallelIO

//------------------------------------------------------------------------------
// Destroy parallel decompositions
void HorzMesh::finalizeParallelIO() {

   IO::destroyDecomp(CellDecompR8);
   IO::destroyDecomp(EdgeDecompR8);
   IO::destroyDecomp(VertexDecompR8);
   IO::destroyDecomp(OnEdgeDecompR8);
   IO::destroyDecomp(OnVertexDecompR8);

} // end finalizeParallelIO

// Read 1D vertex array
void HorzMesh::readVertexArray(HostArray1DReal &VertexArrayH,
                               const std::string &MPASName) {
   Error Err;

   std::string OmegaName;
   std::transform(MPASName.begin(), MPASName.end(), OmegaName.begin(),
                  [](unsigned char c) { return std::toupper(c); });

   // Temporary double precision array for reading
   HostArray1DR8 TmpArrayR8(OmegaName + "Tmp", NVerticesSize);
   int ArrayID;
   Err = IO::readArray(TmpArrayR8.data(), NVerticesAll, MPASName, MeshFileID,
                       VertexDecompR8, ArrayID);
   CHECK_ERROR_ABORT(Err, "HorzMesh: error reading {}", MPASName);

   // Create host array of desired precision and copy the read data into it
   VertexArrayH = HostArray1DReal(OmegaName + "H", NVerticesSize);
   deepCopy(VertexArrayH, TmpArrayR8);
}

// Read 1D edge array
void HorzMesh::readEdgeArray(HostArray1DReal &EdgeArrayH,
                             const std::string &MPASName) {
   Error Err;

   std::string OmegaName;
   std::transform(MPASName.begin(), MPASName.end(), OmegaName.begin(),
                  [](unsigned char c) { return std::toupper(c); });

   // Temporary double precision array for reading
   HostArray1DR8 TmpArrayR8(OmegaName + "Tmp", NEdgesSize);
   int ArrayID;
   Err = IO::readArray(TmpArrayR8.data(), NEdgesAll, MPASName, MeshFileID,
                       EdgeDecompR8, ArrayID);
   CHECK_ERROR_ABORT(Err, "HorzMesh: error reading {}", MPASName);

   // Create host array of desired precision and copy the read data into it
   EdgeArrayH = HostArray1DReal(OmegaName + "H", NEdgesSize);
   deepCopy(EdgeArrayH, TmpArrayR8);
}

// Read 1D cell array
void HorzMesh::readCellArray(HostArray1DReal &CellArrayH,
                             const std::string &MPASName) {
   Error Err;

   std::string OmegaName;
   std::transform(MPASName.begin(), MPASName.end(), OmegaName.begin(),
                  [](unsigned char c) { return std::toupper(c); });

   // Temporary double precision array for reading
   HostArray1DR8 TmpArrayR8(OmegaName + "Tmp", NCellsSize);
   int ArrayID;
   Err = IO::readArray(TmpArrayR8.data(), NCellsAll, MPASName, MeshFileID,
                       CellDecompR8, ArrayID);
   CHECK_ERROR_ABORT(Err, "HorzMesh: error reading {}", MPASName);

   // Create host array of desired precision and copy the read data into it
   CellArrayH = HostArray1DReal(OmegaName + "H", NCellsSize);
   deepCopy(CellArrayH, TmpArrayR8);
}

//------------------------------------------------------------------------------
// Read x/y/z and lon/lat coordinates for cells, edges, and vertices
void HorzMesh::readCoordinates() {

   // Read mesh cell coordinates
   readCellArray(XCellH, "xCell");
   readCellArray(YCellH, "yCell");
   readCellArray(ZCellH, "zCell");

   readCellArray(LonCellH, "lonCell");
   readCellArray(LatCellH, "latCell");

   // Read mesh edge coordinate
   readEdgeArray(XEdgeH, "xEdge");
   readEdgeArray(YEdgeH, "yEdge");
   readEdgeArray(ZEdgeH, "zEdge");

   readEdgeArray(LonEdgeH, "lonEdge");
   readEdgeArray(LatEdgeH, "latEdge");

   // Read mesh vertex coordinates
   readVertexArray(XVertexH, "xVertex");
   readVertexArray(YVertexH, "yVertex");
   readVertexArray(ZVertexH, "zVertex");

   readVertexArray(LonVertexH, "lonVertex");
   readVertexArray(LatVertexH, "latVertex");

} // end readCoordinates

//------------------------------------------------------------------------------
// Read the cell-centered bottom depth
void HorzMesh::readBottomDepth() {
   readCellArray(BottomDepthH, "bottomDepth");
} // end readDepth

//------------------------------------------------------------------------------
// Read the mesh areas (cell, triangle, and kite),
// lengths (between centers and vertices), and edge angles
void HorzMesh::readMeasurements() {

   readCellArray(AreaCellH, "areaCell");

   readVertexArray(AreaTriangleH, "areaTriangle");

   readEdgeArray(DvEdgeH, "dvEdge");

   readEdgeArray(DcEdgeH, "dcEdge");

   readEdgeArray(AngleEdgeH, "angleEdge");

   readCellArray(MeshDensityH, "meshDensity");

   // not using helper function since it kiteAreas is a 2d array
   Error Err;

   // Read into a temporary double precision array
   int KiteAreasOnVertexID;

   HostArray2DR8 TmpKiteAreasOnVertexR8("KiteAreasOnVertex", NVerticesSize,
                                        VertexDegree);
   Err = IO::readArray(TmpKiteAreasOnVertexR8.data(),
                       NVerticesAll * VertexDegree, "kiteAreasOnVertex",
                       MeshFileID, OnVertexDecompR8, KiteAreasOnVertexID);
   CHECK_ERROR_ABORT(Err, "HorzMesh: error reading kiteAreasOnVertex");

   // Create and fill array with Real precision
   KiteAreasOnVertexH =
       HostArray2DReal("KiteAreasOnVertex", NVerticesSize, VertexDegree);
   deepCopy(KiteAreasOnVertexH, TmpKiteAreasOnVertexR8);

} // end readMeasurements

//------------------------------------------------------------------------------
// Read the edge weights used in the discrete potential vorticity flux term
void HorzMesh::readWeights() {

   Error Err;

   int WeightsOnEdgeID;
   HostArray2DR8 TmpWeightsOnEdgeR8("WeightsOnEdge", NEdgesSize, MaxEdges2);
   Err = IO::readArray(TmpWeightsOnEdgeR8.data(), NEdgesAll * MaxEdges2,
                       "weightsOnEdge", MeshFileID, OnEdgeDecompR8,
                       WeightsOnEdgeID);
   CHECK_ERROR_ABORT(Err, "HorzMesh: error reading weightsOnEdge");

   WeightsOnEdgeH = HostArray2DReal("WeightsOnEdge", NEdgesSize, MaxEdges2);
   deepCopy(WeightsOnEdgeH, TmpWeightsOnEdgeR8);

} // end readWeights

//------------------------------------------------------------------------------
// Read the Coriolis parameter at the cells, edges, and vertices
void HorzMesh::readCoriolis() {

   readCellArray(FCellH, "fCell");

   readVertexArray(FVertexH, "fVertex");

   readEdgeArray(FEdgeH, "fEdge");

} // end readCoriolis

//------------------------------------------------------------------------------
// Compute the sign of edge contributions to a cell/vertex for each edge
void HorzMesh::computeEdgeSign() {

   EdgeSignOnCell = Array2DReal("EdgeSignOnCell", NCellsSize, MaxEdges);

   OMEGA_SCOPE(o_NEdgesOnCell, NEdgesOnCell);
   OMEGA_SCOPE(o_EdgesOnCell, EdgesOnCell);
   OMEGA_SCOPE(o_CellsOnEdge, CellsOnEdge);
   OMEGA_SCOPE(o_EdgeSignOnCell, EdgeSignOnCell);

   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int Cell) {
          for (int i = 0; i < o_NEdgesOnCell(Cell); i++) {
             int Edge = o_EdgesOnCell(Cell, i);

             // Vector points from cell 0 to cell 1
             if (Cell == o_CellsOnEdge(Edge, 0)) {
                o_EdgeSignOnCell(Cell, i) = -1.0;
             } else {
                o_EdgeSignOnCell(Cell, i) = 1.0;
             }
          }
       });

   EdgeSignOnCellH = createHostMirrorCopy(EdgeSignOnCell);

   EdgeSignOnVertex =
       Array2DReal("EdgeSignOnVertex", NVerticesSize, VertexDegree);

   OMEGA_SCOPE(o_VertexDegree, VertexDegree);
   OMEGA_SCOPE(o_EdgesOnVertex, EdgesOnVertex);
   OMEGA_SCOPE(o_VerticesOnEdge, VerticesOnEdge);
   OMEGA_SCOPE(o_EdgeSignOnVertex, EdgeSignOnVertex);

   parallelFor(
       {NVerticesAll}, KOKKOS_LAMBDA(int Vertex) {
          for (int i = 0; i < o_VertexDegree; i++) {
             int Edge = o_EdgesOnVertex(Vertex, i);

             // Vector points from vertex 0 to vertex 1
             if (Vertex == o_VerticesOnEdge(Edge, 0)) {
                o_EdgeSignOnVertex(Vertex, i) = -1.0;
             } else {
                o_EdgeSignOnVertex(Vertex, i) = 1.0;
             }
          }
       });

   EdgeSignOnVertexH = createHostMirrorCopy(EdgeSignOnVertex);
} // end computeEdgeSign

//------------------------------------------------------------------------------
// Set mesh scaling coefficients for mixing terms in momentum and tracer
// equations so viscosity and diffusion scale with mesh.
void HorzMesh::setMeshScaling() {

   MeshScalingDel2 = Array1DReal("MeshScalingDel2", NEdgesSize);
   MeshScalingDel4 = Array1DReal("MeshScalingDel4", NEdgesSize);

   OMEGA_SCOPE(o_MeshScalingDel2, MeshScalingDel2);
   OMEGA_SCOPE(o_MeshScalingDel4, MeshScalingDel4);

   // TODO: implement mesh scaling by cell area, only no scaling
   // option for now
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int Edge) {
          o_MeshScalingDel2(Edge) = 1.0;
          o_MeshScalingDel4(Edge) = 1.0;
       });

   MeshScalingDel2H = createHostMirrorCopy(MeshScalingDel2);
   MeshScalingDel4H = createHostMirrorCopy(MeshScalingDel4);

} // end setMeshScaling

//------------------------------------------------------------------------------
// Perform copy to device for mesh variables
void HorzMesh::copyToDevice() {

   AreaCell          = createDeviceMirrorCopy(AreaCellH);
   AreaTriangle      = createDeviceMirrorCopy(AreaTriangleH);
   KiteAreasOnVertex = createDeviceMirrorCopy(KiteAreasOnVertexH);
   DcEdge            = createDeviceMirrorCopy(DcEdgeH);
   DvEdge            = createDeviceMirrorCopy(DvEdgeH);
   AngleEdge         = createDeviceMirrorCopy(AngleEdgeH);
   WeightsOnEdge     = createDeviceMirrorCopy(WeightsOnEdgeH);
   FVertex           = createDeviceMirrorCopy(FVertexH);
   BottomDepth       = createDeviceMirrorCopy(BottomDepthH);
   FEdge             = createDeviceMirrorCopy(FEdgeH);
   XCell             = createDeviceMirrorCopy(XCellH);
   YCell             = createDeviceMirrorCopy(YCellH);
   XEdge             = createDeviceMirrorCopy(XEdgeH);
   YEdge             = createDeviceMirrorCopy(YEdgeH);

} // end copyToDevice

//------------------------------------------------------------------------------
// Get default mesh
HorzMesh *HorzMesh::getDefault() { return HorzMesh::DefaultHorzMesh; }

//------------------------------------------------------------------------------
// Get mesh by name
HorzMesh *HorzMesh::get(const std::string Name ///< [in] Name of mesh
) {

   // look for an instance of this name
   auto it = AllHorzMeshes.find(Name);

   // if found, return the mesh pointer
   if (it != AllHorzMeshes.end()) {
      return it->second.get();

      // otherwise print error and return null pointer
   } else {
      LOG_ERROR("HorzMesh::get: Attempt to retrieve non-existent mesh:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }
} // end get mesh

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
