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
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"

namespace OMEGA {

// create the static class members
HorzMesh *HorzMesh::DefaultHorzMesh = nullptr;
std::map<std::string, std::unique_ptr<HorzMesh>> HorzMesh::AllHorzMeshes;

//------------------------------------------------------------------------------
// Initialize the mesh. Assumes that Decomp has already been initialized.

int HorzMesh::init() {

   int Err = 0; // default successful return code

   // Retrieve the default decomposition
   Decomp *DefDecomp = Decomp::getDefault();

   // Create the default mesh and set pointer to it
   HorzMesh::DefaultHorzMesh = create("Default", DefDecomp);
   return Err;
}

//------------------------------------------------------------------------------
// Construct a new local mesh given a decomposition

HorzMesh::HorzMesh(const std::string &Name, //< [in] Name for new mesh
                   Decomp *MeshDecomp       //< [in] Decomp for the new mesh
) {

   MeshName = Name;

   // Retrieve mesh files name from Decomp
   MeshFileName = MeshDecomp->MeshFileName;

   // Retrieve mesh cell/edge/vertex totals from Decomp
   NCellsHalo  = MeshDecomp->NCellsHalo;
   NCellsHaloH = MeshDecomp->NCellsHaloH;
   NCellsOwned = MeshDecomp->NCellsOwned;
   NCellsAll   = MeshDecomp->NCellsAll;
   NCellsSize  = MeshDecomp->NCellsSize;

   NEdgesHalo     = MeshDecomp->NEdgesHalo;
   NEdgesHaloH    = MeshDecomp->NEdgesHaloH;
   NEdgesOwned    = MeshDecomp->NEdgesOwned;
   NEdgesAll      = MeshDecomp->NEdgesAll;
   NEdgesSize     = MeshDecomp->NEdgesSize;
   MaxCellsOnEdge = MeshDecomp->MaxCellsOnEdge;
   MaxEdges       = MeshDecomp->MaxEdges;

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
   I4 Err;
   Err = OMEGA::IO::openFile(MeshFileID, MeshFileName, IO::ModeRead);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error opening mesh file");

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

   // TODO: implement setMasks during Mesh constructor
   // setMasks()

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

   I4 Err;
   I4 NDims             = 1;
   IO::Rearranger Rearr = IO::RearrBox;

   // Create the IO decomp for arrays with (NCells) dimensions
   std::vector<I4> CellDims{MeshDecomp->NCellsGlobal};
   std::vector<I4> CellID(NCellsAll);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      CellID[Cell] = MeshDecomp->CellIDH(Cell) - 1;
   }

   Err = IO::createDecomp(CellDecompR8, IO::IOTypeR8, NDims, CellDims,
                          NCellsAll, CellID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating cell IO decomposition");

   // Create the IO decomp for arrays with (NEdges) dimensions
   std::vector<I4> EdgeDims{MeshDecomp->NEdgesGlobal};
   std::vector<I4> EdgeID(NEdgesAll);
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      EdgeID[Edge] = MeshDecomp->EdgeIDH(Edge) - 1;
   }

   Err = IO::createDecomp(EdgeDecompR8, IO::IOTypeR8, NDims, EdgeDims,
                          NEdgesAll, EdgeID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating edge IO decomposition");

   // Create the IO decomp for arrays with (NVertices) dimensions
   std::vector<I4> VertexDims{MeshDecomp->NVerticesGlobal};
   std::vector<I4> VertexID(NVerticesAll);
   for (int Vertex = 0; Vertex < NVerticesAll; ++Vertex) {
      VertexID[Vertex] = MeshDecomp->VertexIDH(Vertex) - 1;
   }

   Err = IO::createDecomp(VertexDecompR8, IO::IOTypeR8, NDims, VertexDims,
                          NVerticesAll, VertexID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating vertex IO decomposition");

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

   Err = IO::createDecomp(OnEdgeDecompR8, IO::IOTypeR8, NDims, OnEdgeDims2,
                          OnEdgeSize2, OnEdgeOffset2, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating OnEdge IO decomposition");

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

   Err = IO::createDecomp(OnVertexDecompR8, IO::IOTypeR8, NDims, OnVertexDims,
                          OnVertexSize, OnVertexOffset, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating OnVertex IO decomposition");

} // end initParallelIO

//------------------------------------------------------------------------------
// Destroy parallel decompositions
void HorzMesh::finalizeParallelIO() {

   int Err = 0; // default return code

   // Destroy the IO decomp for arrays with (NCells) dimensions
   Err = IO::destroyDecomp(CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error destroying cell IO decomposition");

   // Destroy the IO decomp for arrays with (NEdges) dimensions
   Err = IO::destroyDecomp(EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error destroying edge IO decomposition");

   // Destroy the IO decomp for arrays with (NVertices) dimensions
   Err = IO::destroyDecomp(VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error destroying vertex IO decomposition");

   // Destroy the IO decomp for arrays with (NEdges, 2*MaxEdges) dimensions
   Err = IO::destroyDecomp(OnEdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error destroying OnEdge IO decomposition");

   // Destroy the IO decomp for arrays with (NVertices, VertexDegree) dimensions
   Err = IO::destroyDecomp(OnVertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error destroying OnVertex IO decomposition");

} // end finalizeParallelIO

//------------------------------------------------------------------------------
// Read x/y/z and lon/lat coordinates for cells, edges, and vertices
void HorzMesh::readCoordinates() {

   I4 Err;

   // Read mesh cell coordinates
   int XCellID;
   XCellH = HostArray1DR8("XCell", NCellsSize);
   Err    = IO::readArray(XCellH.data(), NCellsAll, "xCell", MeshFileID,
                          CellDecompR8, XCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xCell");

   int YCellID;
   YCellH = HostArray1DR8("YCell", NCellsSize);
   Err    = IO::readArray(YCellH.data(), NCellsAll, "yCell", MeshFileID,
                          CellDecompR8, YCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yCell");

   int ZCellID;
   ZCellH = HostArray1DR8("ZCell", NCellsSize);
   Err    = IO::readArray(ZCellH.data(), NCellsAll, "zCell", MeshFileID,
                          CellDecompR8, ZCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zCell");

   int LonCellID;
   LonCellH = HostArray1DR8("LonCell", NCellsSize);
   Err      = IO::readArray(LonCellH.data(), NCellsAll, "lonCell", MeshFileID,
                            CellDecompR8, LonCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonCell");

   int LatCellID;
   LatCellH = HostArray1DR8("LatCell", NCellsSize);
   Err      = IO::readArray(LatCellH.data(), NCellsAll, "latCell", MeshFileID,
                            CellDecompR8, LatCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latCell");

   // Read mesh edge coordinateID
   int XEdgeID;
   XEdgeH = HostArray1DR8("XEdge", NEdgesSize);
   Err    = IO::readArray(XEdgeH.data(), NEdgesAll, "xEdge", MeshFileID,
                          EdgeDecompR8, XEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xEdge");

   int YEdgeID;
   YEdgeH = HostArray1DR8("YEdge", NEdgesSize);
   Err    = IO::readArray(YEdgeH.data(), NEdgesAll, "yEdge", MeshFileID,
                          EdgeDecompR8, YEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yEdge");

   int ZEdgeID;
   ZEdgeH = HostArray1DR8("ZEdge", NEdgesSize);
   Err    = IO::readArray(ZEdgeH.data(), NEdgesAll, "zEdge", MeshFileID,
                          EdgeDecompR8, ZEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zEdge");

   int LonEdgeID;
   LonEdgeH = HostArray1DR8("LonEdge", NEdgesSize);
   Err      = IO::readArray(LonEdgeH.data(), NEdgesAll, "lonEdge", MeshFileID,
                            EdgeDecompR8, LonEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonEdge");

   int LatEdgeID;
   LatEdgeH = HostArray1DR8("LatEdge", NEdgesSize);
   Err      = IO::readArray(LatEdgeH.data(), NEdgesAll, "latEdge", MeshFileID,
                            EdgeDecompR8, LatEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latEdge");

   // Read mesh vertex coordinates
   int XVertexID;
   XVertexH = HostArray1DR8("XVertex", NVerticesSize);
   Err = IO::readArray(XVertexH.data(), NVerticesAll, "xVertex", MeshFileID,
                       VertexDecompR8, XVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xVertex");

   int YVertexID;
   YVertexH = HostArray1DR8("YVertex", NVerticesSize);
   Err = IO::readArray(YVertexH.data(), NVerticesAll, "yVertex", MeshFileID,
                       VertexDecompR8, YVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yVertex");

   int ZVertexID;
   ZVertexH = HostArray1DR8("ZVertex", NVerticesSize);
   Err = IO::readArray(ZVertexH.data(), NVerticesAll, "zVertex", MeshFileID,
                       VertexDecompR8, ZVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zVertex");

   int LonVertexID;
   LonVertexH = HostArray1DR8("LonVertex", NVerticesSize);
   Err = IO::readArray(LonVertexH.data(), NVerticesAll, "lonVertex", MeshFileID,
                       VertexDecompR8, LonVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonVertex");

   int LatVertexID;
   LatVertexH = HostArray1DR8("LatVertex", NVerticesSize);
   Err = IO::readArray(LatVertexH.data(), NVerticesAll, "latVertex", MeshFileID,
                       VertexDecompR8, LatVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latVertex");

} // end readCoordinates

//------------------------------------------------------------------------------
// Read the cell-centered bottom depth
void HorzMesh::readBottomDepth() {

   I4 Err;

   int BottomDepthID;
   BottomDepthH = HostArray1DR8("BottomDepth", NCellsSize);
   Err          = IO::readArray(BottomDepthH.data(), NCellsAll, "bottomDepth",
                                MeshFileID, CellDecompR8, BottomDepthID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading bottomDepth");

} // end readDepth

//------------------------------------------------------------------------------
// Read the mesh areas (cell, triangle, and kite),
// lengths (between centers and vertices), and edge angles
void HorzMesh::readMeasurements() {

   I4 Err;

   int AreaCellID;
   AreaCellH = HostArray1DR8("AreaCell", NCellsSize);
   Err = IO::readArray(AreaCellH.data(), NCellsAll, "areaCell", MeshFileID,
                       CellDecompR8, AreaCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading areaCell");

   int AreaTriangleID;
   AreaTriangleH = HostArray1DR8("AreaTriangle", NVerticesSize);
   Err = IO::readArray(AreaTriangleH.data(), NVerticesAll, "areaTriangle",
                       MeshFileID, VertexDecompR8, AreaTriangleID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading areaTriangle");

   int DvEdgeID;
   DvEdgeH = HostArray1DR8("DvEdge", NEdgesSize);
   Err     = IO::readArray(DvEdgeH.data(), NEdgesAll, "dvEdge", MeshFileID,
                           EdgeDecompR8, DvEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading dvEdge");

   int DcEdgeID;
   DcEdgeH = HostArray1DR8("DcEdge", NEdgesSize);
   Err     = IO::readArray(DcEdgeH.data(), NEdgesAll, "dcEdge", MeshFileID,
                           EdgeDecompR8, DcEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading dcEdge");

   int AngleEdgeID;
   AngleEdgeH = HostArray1DR8("AngleEdge", NEdgesSize);
   Err = IO::readArray(AngleEdgeH.data(), NEdgesAll, "angleEdge", MeshFileID,
                       EdgeDecompR8, AngleEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading angleEdge");

   int MeshDensityID;
   MeshDensityH = HostArray1DR8("MeshDensity", NCellsSize);
   Err          = IO::readArray(MeshDensityH.data(), NCellsAll, "meshDensity",
                                MeshFileID, CellDecompR8, MeshDensityID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading meshDensity");

   int KiteAreasOnVertexID;
   KiteAreasOnVertexH =
       HostArray2DR8("KiteAreasOnVertex", NVerticesSize, VertexDegree);
   Err = IO::readArray(KiteAreasOnVertexH.data(), NVerticesAll * VertexDegree,
                       "kiteAreasOnVertex", MeshFileID, OnVertexDecompR8,
                       KiteAreasOnVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading kiteAreasOnVertex");

} // end readMeasurements

//------------------------------------------------------------------------------
// Read the edge weights used in the discrete potential vorticity flux term
void HorzMesh::readWeights() {

   I4 Err;

   int WeightsOnEdgeID;
   WeightsOnEdgeH = HostArray2DR8("WeightsOnEdge", NEdgesSize, MaxEdges2);
   Err            = IO::readArray(WeightsOnEdgeH.data(), NEdgesAll * MaxEdges2,
                                  "weightsOnEdge", MeshFileID, OnEdgeDecompR8,
                                  WeightsOnEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading weightsOnEdge");

} // end readWeights

//------------------------------------------------------------------------------
// Read the Coriolis parameter at the cells, edges, and vertices
void HorzMesh::readCoriolis() {

   int Err;

   int FCellID;
   FCellH = HostArray1DR8("FCell", NCellsSize);
   Err    = IO::readArray(FCellH.data(), NCellsAll, "fCell", MeshFileID,
                          CellDecompR8, FCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading fCell");

   int FVertexID;
   FVertexH = HostArray1DR8("FVertex", NVerticesSize);
   Err = IO::readArray(FVertexH.data(), NVerticesAll, "fVertex", MeshFileID,
                       VertexDecompR8, FVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading fVertex");

   int FEdgeID;
   FEdgeH = HostArray1DR8("FEdge", NEdgesSize);
   Err    = IO::readArray(FEdgeH.data(), NEdgesAll, "fEdge", MeshFileID,
                          EdgeDecompR8, FEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading fEdge");

} // end readCoriolis

//------------------------------------------------------------------------------
// Compute the sign of edge contributions to a cell/vertex for each edge
void HorzMesh::computeEdgeSign() {

   EdgeSignOnCell = Array2DR8("EdgeSignOnCell", NCellsSize, MaxEdges);

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
       Array2DR8("EdgeSignOnVertex", NVerticesSize, VertexDegree);

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
// set computational masks for mesh elements
// TODO: this is just a placeholder, implement actual masks for edges, cells,
// and vertices
void HorzMesh::setMasks(int NVertLevels) {

   EdgeMask = Array2DR8("EdgeMask", NEdgesSize, NVertLevels);

   OMEGA_SCOPE(O_EdgeMask, EdgeMask);

   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int Edge) {
          for (int K = 0; K < NVertLevels; ++K) {
             O_EdgeMask(Edge, K) = 1.0;
          }
       });

   EdgeMaskH = createHostMirrorCopy(EdgeMask);

} // end setMasks

//------------------------------------------------------------------------------
// Set mesh scaling coefficients for mixing terms in momentum and tracer
// equations so viscosity and diffusion scale with mesh.
void HorzMesh::setMeshScaling() {

   MeshScalingDel2 = Array1DR8("MeshScalingDel2", NEdgesSize);
   MeshScalingDel4 = Array1DR8("MeshScalingDel4", NEdgesSize);

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
