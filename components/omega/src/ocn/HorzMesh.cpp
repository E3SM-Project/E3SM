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
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"

namespace OMEGA {
//------------------------------------------------------------------------------
// Initialize the mesh 

int HorzMesh::init() { 

   int Err = 0; // default successful return code

   // Retrieve the default decomposition
   Decomp *DefDecomp = Decomp::getDefault();

   // Create the default mesh 
   HorzMesh DefHorzMesh("Default", DefDecomp);

   // Retrieve this mesh and set pointer to DefaultHorzMesh
   HorzMesh::DefaultHorzMesh = HorzMesh::get("Default");

   return Err;
}

//------------------------------------------------------------------------------
// Construct a new local mesh given a decomposition

HorzMesh::HorzMesh(const std::string & Name, //< [in] Name for new mesh
                   Decomp *MeshDecomp        //< [in] Decomp for the new mesh
) {

   // Retrieve mesh files name from Decomp
   MeshFileName = MeshDecomp->MeshFileName;

   // Retrieve mesh cell/edge/vertex totals from Decomp
   NCellsOwned = MeshDecomp->NCellsOwned;
   NCellsAll   = MeshDecomp->NCellsAll;
   NCellsSize  = MeshDecomp->NCellsSize;

   NEdgesOwned    = MeshDecomp->NEdgesOwned;
   NEdgesAll      = MeshDecomp->NEdgesAll;
   NEdgesSize     = MeshDecomp->NEdgesSize;
   MaxCellsOnEdge = MeshDecomp->MaxCellsOnEdge;
   MaxEdges       = MeshDecomp->MaxEdges;

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

   // Copy host data to device
   copyToDevice();

   // TODO: add ability to compute (rather than read in)
   // dependent mesh quantities

   // Compute EdgeSignOnCells and EdgeSignOnVertex
   computeEdgeSign();

   // Associate this instance with a name 
   AllHorzMeshes.emplace(Name, *this);

} // end horizontal mesh constructor

//------------------------------------------------------------------------------
// Destroys a local mesh and deallocates all arrays
HorzMesh::~HorzMesh() {

   // TODO: add deletes for all arrays and remove from AllDecomps map

} // end deconstructor

//------------------------------------------------------------------------------
// Deallocates arrays
void HorzMesh::clear() {

   AreaCell.deallocate();
   AreaTriangle.deallocate();
   KiteAreasOnVertex.deallocate();
   DvEdge.deallocate();
   DcEdge.deallocate();
   AngleEdge.deallocate();
   WeightsOnEdge.deallocate();
   FVertex.deallocate();
   BottomDepth.deallocate();
   EdgeSignOnCell.deallocate();
   EdgeSignOnCellH.deallocate();
   EdgeSignOnVertex.deallocate();
   EdgeSignOnVertexH.deallocate();

} // end clear

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
 

//------------------------------------------------------------------------------
// Read x/y/z and lon/lat coordinates for cells, edges, and vertices
void HorzMesh::readCoordinates() {

   I4 Err;

   // Read mesh cell coordinates
   int XCellID;
   XCellH = ArrayHost1DR8("XCell", NCellsAll);
   Err    = IO::readArray(XCellH.data(), NCellsAll, "xCell", MeshFileID,
                          CellDecompR8, XCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xCell");

   int YCellID;
   YCellH = ArrayHost1DR8("YCell", NCellsAll);
   Err    = IO::readArray(YCellH.data(), NCellsAll, "yCell", MeshFileID,
                          CellDecompR8, YCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yCell");

   int ZCellID;
   ZCellH = ArrayHost1DR8("ZCell", NCellsAll);
   Err    = IO::readArray(ZCellH.data(), NCellsAll, "zCell", MeshFileID,
                          CellDecompR8, ZCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zCell");

   int LonCellID;
   LonCellH = ArrayHost1DR8("LonCell", NCellsAll);
   Err      = IO::readArray(LonCellH.data(), NCellsAll, "lonCell", MeshFileID,
                            CellDecompR8, LonCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonCell");

   int LatCellID;
   LatCellH = ArrayHost1DR8("LatCell", NCellsAll);
   Err      = IO::readArray(LatCellH.data(), NCellsAll, "latCell", MeshFileID,
                            CellDecompR8, LatCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latCell");

   // Read mesh edge coordinateID
   int XEdgeID;
   XEdgeH = ArrayHost1DR8("XEdge", NEdgesAll);
   Err    = IO::readArray(XEdgeH.data(), NEdgesAll, "xEdge", MeshFileID,
                          EdgeDecompR8, XEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xEdge");

   int YEdgeID;
   YEdgeH = ArrayHost1DR8("YEdge", NEdgesAll);
   Err    = IO::readArray(YEdgeH.data(), NEdgesAll, "yEdge", MeshFileID,
                          EdgeDecompR8, YEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yEdge");

   int ZEdgeID;
   ZEdgeH = ArrayHost1DR8("ZEdge", NEdgesAll);
   Err    = IO::readArray(ZEdgeH.data(), NEdgesAll, "zEdge", MeshFileID,
                          EdgeDecompR8, ZEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zEdge");

   int LonEdgeID;
   LonEdgeH = ArrayHost1DR8("LonEdge", NEdgesAll);
   Err      = IO::readArray(LonEdgeH.data(), NEdgesAll, "lonEdge", MeshFileID,
                            EdgeDecompR8, LonEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonEdge");

   int LatEdgeID;
   LatEdgeH = ArrayHost1DR8("LatEdge", NEdgesAll);
   Err      = IO::readArray(LatEdgeH.data(), NEdgesAll, "latEdge", MeshFileID,
                            EdgeDecompR8, LatEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latEdge");

   // Read mesh vertex coordinates
   int XVertexID;
   XVertexH = ArrayHost1DR8("XVertex", NVerticesAll);
   Err = IO::readArray(XVertexH.data(), NVerticesAll, "xVertex", MeshFileID,
                       VertexDecompR8, XVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xVertex");

   int YVertexID;
   YVertexH = ArrayHost1DR8("YVertex", NVerticesAll);
   Err = IO::readArray(YVertexH.data(), NVerticesAll, "yVertex", MeshFileID,
                       VertexDecompR8, YVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yVertex");

   int ZVertexID;
   ZVertexH = ArrayHost1DR8("ZVertex", NVerticesAll);
   Err = IO::readArray(ZVertexH.data(), NVerticesAll, "zVertex", MeshFileID,
                       VertexDecompR8, ZVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zVertex");

   int LonVertexID;
   LonVertexH = ArrayHost1DR8("LonVertex", NVerticesAll);
   Err = IO::readArray(LonVertexH.data(), NVerticesAll, "lonVertex", MeshFileID,
                       VertexDecompR8, LonVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonVertex");

   int LatVertexID;
   LatVertexH = ArrayHost1DR8("LatVertex", NVerticesAll);
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
   BottomDepthH = ArrayHost1DR8("BottomDepth", NCellsAll);
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
   AreaCellH = ArrayHost1DR8("AreaCell", NCellsAll);
   Err = IO::readArray(AreaCellH.data(), NCellsAll, "areaCell", MeshFileID,
                       CellDecompR8, AreaCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading areaCell");

   int AreaTriangleID;
   AreaTriangleH = ArrayHost1DR8("AreaTriangle", NVerticesAll);
   Err = IO::readArray(AreaTriangleH.data(), NVerticesAll, "areaTriangle",
                       MeshFileID, VertexDecompR8, AreaTriangleID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading areaTriangle");

   int DvEdgeID;
   DvEdgeH = ArrayHost1DR8("DvEdge", NEdgesAll);
   Err     = IO::readArray(DvEdgeH.data(), NEdgesAll, "dvEdge", MeshFileID,
                           EdgeDecompR8, DvEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading dvEdge");

   int DcEdgeID;
   DcEdgeH = ArrayHost1DR8("DcEdge", NEdgesAll);
   Err     = IO::readArray(DcEdgeH.data(), NEdgesAll, "dcEdge", MeshFileID,
                           EdgeDecompR8, DcEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading dcEdge");

   int AngleEdgeID;
   AngleEdgeH = ArrayHost1DR8("AngleEdge", NEdgesAll);
   Err = IO::readArray(AngleEdgeH.data(), NEdgesAll, "angleEdge", MeshFileID,
                       EdgeDecompR8, AngleEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading angleEdge");

   int MeshDensityID;
   MeshDensityH = ArrayHost1DR8("MeshDensity", NCellsAll);
   Err          = IO::readArray(MeshDensityH.data(), NCellsAll, "meshDensity",
                                MeshFileID, CellDecompR8, MeshDensityID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading meshDensity");

   int KiteAreasOnVertexID;
   KiteAreasOnVertexH =
       ArrayHost2DR8("KiteAreasOnVertex", NVerticesAll, VertexDegree);
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
   WeightsOnEdgeH = ArrayHost2DR8("WeightsOnEdge", NEdgesAll, MaxEdges2);
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
   FCellH = ArrayHost1DR8("FCell", NCellsAll);
   Err    = IO::readArray(FCellH.data(), NCellsAll, "fCell", MeshFileID,
                          CellDecompR8, FCellID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading fCell");

   int FVertexID;
   FVertexH = ArrayHost1DR8("FVertex", NVerticesAll);
   Err = IO::readArray(FVertexH.data(), NVerticesAll, "fVertex", MeshFileID,
                       VertexDecompR8, FVertexID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading fVertex");

   int FEdgeID;
   FEdgeH = ArrayHost1DR8("FEdge", NEdgesAll);
   Err    = IO::readArray(FEdgeH.data(), NEdgesAll, "fEdge", MeshFileID,
                          EdgeDecompR8, FEdgeID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading fEdge");

} // end readCoriolis

//------------------------------------------------------------------------------
// Compute the sign of edge contributions to a cell/vertex for each edge
void HorzMesh::computeEdgeSign() {

   EdgeSignOnCell = Array2DR8("EdgeSignOnCell", NCellsAll, MaxEdges);
   yakl::c::parallel_for(
       yakl::c::SimpleBounds<1>(NCellsAll), YAKL_LAMBDA(int Cell) {
          for (int i = 0; i < NEdgesOnCell(Cell); i++) {
             int Edge = EdgesOnCell(Cell, i);

             // Vector points from cell 0 to cell 1
             if (Cell == CellsOnEdge(Edge, 0)) {
                EdgeSignOnCell(Cell, i) = -1.0;
             } else {
                EdgeSignOnCell(Cell, i) = 1.0;
             }
          }
       });
   EdgeSignOnCellH = EdgeSignOnCell.createHostCopy();

   EdgeSignOnVertex = Array2DR8("EdgeSignOnVertex", NVerticesAll, VertexDegree);
   yakl::c::parallel_for(
       yakl::c::SimpleBounds<1>(NVerticesAll), YAKL_LAMBDA(int Vertex) {
          for (int i = 0; i < VertexDegree; i++) {
             int Edge = EdgesOnVertex(Vertex, i);

             // Vector points from vertex 0 to vertex 1
             if (Vertex == VerticesOnEdge(Edge, 0)) {
                EdgeSignOnVertex(Vertex, i) = -1.0;
             } else {
                EdgeSignOnVertex(Vertex, i) = 1.0;
             }
          }
       });
   EdgeSignOnVertexH = EdgeSignOnVertex.createHostCopy();

} // end computeEdgeSign

//------------------------------------------------------------------------------
// Perform copy to device for mesh variables
void HorzMesh::copyToDevice() {

   AreaCell          = AreaCellH.createDeviceCopy();
   AreaTriangle      = AreaTriangleH.createDeviceCopy();
   KiteAreasOnVertex = KiteAreasOnVertexH.createDeviceCopy();
   DcEdge            = DcEdgeH.createDeviceCopy();
   DvEdge            = DvEdgeH.createDeviceCopy();
   AngleEdge         = AngleEdgeH.createDeviceCopy();
   WeightsOnEdge     = WeightsOnEdgeH.createDeviceCopy();
   FVertex           = FVertexH.createDeviceCopy();
   BottomDepth       = BottomDepthH.createDeviceCopy();

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
      return &(it->second);

   // otherwise print error and return null pointer
   } else {
      LOG_ERROR("HorzMesh::get: Attempt to retrieve non-existent mesh:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }
} // end get mesh

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
